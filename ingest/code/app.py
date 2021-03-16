from Bio import SeqIO
from rdflib import BNode, ConjunctiveGraph, Graph, Literal, Namespace, RDF, URIRef, plugin
from rdflib.namespace import XSD
from rdflib.serializer import Serializer
import datapackage
import tableschema.exceptions
from tableschema import Table
from datetime import timedelta
import uuid
import yaml
from minio import Minio
from minio.error import S3Error

def getOntologyMappings(config_file):
    """Get the mappings from the data templates to the ontology"""
    with open(config_file, 'r') as yamlfile:
        mappings = yaml.load(yamlfile, Loader=yaml.FullLoader)
    return mappings

def processFieldValue( value, descriptor, field_type):
    """Process a field's value"""

    if 'missingValues' in descriptor:
        for miss in descriptor['missingValues']:
            if value == miss:
              return None

    #Convert to correct ES data type
    if descriptor['type'] == 'number':
        return float(value)
    elif descriptor['type'] == 'integer':
        return int(value)
    return value

def processField(value, descriptor, field_type, delimiterField='opp:fieldValueDelimiter'):
    """Process a field"""
    if None is value:
        return value

    is_array_value = False
    values = [value]

    # Is the value an array of values?
    if delimiterField in descriptor:
        is_array_value = True
        values = value.split(descriptor[delimiterField])

    # Are there contraints that help process the data?
    if 'constraints' in descriptor:
        if 'pattern' in descriptor['constraints']:
            for idx, val in enumerate(values):
                regex = re.compile('^{0}$'.format(descriptor['constraints']['pattern']))
                match = regex.match(value)
                if match and match.lastindex is 1:
                    #### could handle VALUE PROCESSING here
                    values[idx] = processFieldValue(match.group(1), descriptor, field_type)
                else:
                    # must be null
                    values[idx] = None
    else:
        for idx,val in enumerate(values):
            values[idx] = processFieldValue(val, descriptor, field_type)

    # Return the value
    if is_array_value:
        return values
    else:
        return values[0]

def main():

    with open('./config.yaml', 'r') as yamlfile:
        CONFIG = yaml.load(yamlfile, Loader=yaml.FullLoader)


    # Get the ontology mappings for interpreting the data files
    ontology_mappings = getOntologyMappings('./data-template-ontology-mappings.yaml')

    # Setup the object store client
    mc = Minio(
        CONFIG['s3']['endpoint'],
        access_key=CONFIG['s3']['access_key'],
        secret_key=CONFIG['s3']['secret_key'],
    )

    # Get the datapackage.json for a dataset in the object store
    dp_url = mc.presigned_get_object('ingested', 'Metzyme/data/datapackage.json', expires=timedelta(minutes=1))
    dp = datapackage.DataPackage(dp_url, None)
    if (dp.errors):
        for error in dp.errors:
            logging.error(error)
        raise Exception('Invalid data package')

    # Validate the Datapackage
    try:
        valid = datapackage.validate(dp.descriptor)
    except exceptions.ValidationError as exception:
        for error in datapackage.exception.errors:
            logging.error(error)
        raise Exception('Invalid data package')

    # Define a top-level graph
    cg = ConjunctiveGraph()

    # Define Namespaces
    geo = Namespace("http://www.w3.org/2003/01/geo/wgs84_pos#")
    opp = Namespace("http://schema.oceanproteinportal.org/v2/")
    rdfs = Namespace("http://www.w3.org/2000/01/rdf-schema#")
    schema = Namespace("http://schema.org/")
    sosa = Namespace("http://www.w3.org/ns/sosa/")
    ssn = Namespace("http://www.w3.org/ns/ssn/")
    ssn_system = Namespace("http://www.w3.org/ns/ssn/systems/")
    wikidata = Namespace("http://www.wikidata.org/entity/")

    # Add namespaces to graph
    cg.bind("geo", geo)
    cg.bind("opp", opp)
    cg.bind("rdf", RDF)
    cg.bind("rdfs", rdfs)
    cg.bind("schema", schema)
    cg.bind("sosa", sosa)
    cg.bind("ssn", ssn)
    cg.bind("ssn-system", ssn_system)
    cg.bind("wikidata", wikidata)
    cg.bind("xsd", XSD)

    # Define URI space for OPP subjects
    id_uri_prefix = "http://id.oceanproteinportal.org/"
    odo_dt = "http://ocean-data.org/schema/data-type/"

    # Generate datasetId & define a Dataset resource
    datasetId = str(uuid.uuid3(uuid.NAMESPACE_DNS, dp.descriptor['name'] + '_ver.' + dp.descriptor.get('version', 'noversion')))
    dataset = URIRef(id_uri_prefix + "dataset/" + datasetId)

    # Define a Named Graph
    g = Graph(cg.store, dataset)

    # Construct the Dataset
    g.add((dataset, RDF.type, schema.Dataset))
    # Dataset Identifier
    datasetIdentifier = BNode()
    g.add((datasetIdentifier, RDF.type, opp.Identifier))
    datasetIdScheme = BNode()
    g.add((datasetIdScheme, RDF.type, opp.DatasetIdentfierScheme))
    g.add((datasetIdentifier, opp.identifierScheme, datasetIdScheme))
    g.add((datasetIdentifier, opp.identifierValue, Literal(datasetId, datatype=XSD.token)))
    g.add((dataset, opp.identifier, datasetIdentifier))
    g.add((dataset, opp.definesIdentifierScheme, datasetIdScheme))

    g.add((dataset, schema.name, Literal(dp.descriptor['title'], datatype=XSD.string)))
    g.add((dataset, schema.description, Literal(dp.descriptor['description'], datatype=XSD.string)))
    g.add((dataset, schema.url, Literal(dp.descriptor['homepage'], datatype=XSD.anyURI)))
    g.add((dataset, schema.version, Literal(dp.descriptor['version'], datatype=XSD.token)))
    g.add((dataset, schema.alternateName, Literal(dp.descriptor['opp:shortName'], datatype=XSD.string)))

    for keyword in dp.descriptor['keywords']:
      g.add((dataset, schema.keyword, Literal(keyword, datatype=XSD.token)))

    #TODO: Figure out how to handle the License, if we want to
    #TODO: Add bibliographic citation (dcterms:bibliographicCitation) to DP
    #TODO: Add any other products this dataset cites (schema:citation) to DP

    for index, contributor in enumerate(dp.descriptor['contributors']):
        # Agent
        agent = BNode("agent-{idx:.2f}".format(idx=index))
        if 'uri' in contributor:
            agent = URIRef(contributor['uri'])
            # TODO: SHACL rule to infer proper Agent subclass

        g.add((agent, RDF.type, opp.Agent))
        g.add((agent, schema.name, Literal(contributor['title'], datatype=XSD.string)))
        if 'email' in contributor:
            g.add((agent, schema.email, Literal(contributor['email'], datatype=XSD.string)))

        # Deal with an ORCID
        if 'orcid' in contributor:
            orcid = BNode()
            g.add((orcid, RDF.type, opp.Identifier))
            g.add((orcid, opp.identifierScheme, wikidata.Q51044))
            g.add((orcid, opp.identifierValue, Literal(contributor['orcid'], datatype=XSD.token)))
            g.add((agent, opp.identifier, orcid))

        # Role
        if 'role' in contributor:
            role = opp.Role_Contributor
        if 'author' == contributor['role']:
            role = opp.Role_Author
        elif 'contact' == contributor['role']:
            role = opp.Role_Contact
        elif 'publisher' == contributor['role']:
            role = opp.Role_Publisher

        # AgentRole
        agent_role = BNode()
        g.add((agent_role, RDF.type, opp.AgentRole))
        g.add((agent_role, opp.isAgentRoleFor, dataset))
        g.add((agent_role, opp.performedBy, agent))
        g.add((agent_role, opp.inRole, role))
        # TODO: SHACL rule to infer opp.hasAgentRole

    '''
    TODO - At this point we really don't need these triples in memory anymore.
    1. Can we store them somewhere to be used later and then clear the graph?
    '''
    # Store the dataset RDF
    cg.serialize(destination='./dataset.nq', format='nquads')
    # Clear the in-memory graph
    g.remove((None, None, None))

    #Store any cruises
    cruises = {}
    if 'odo:hasDeployment' in dp.descriptor:
        for deployment in dp.descriptor['odo:hasDeployment']:
            cruises[deployment['name']] = URIRef(deployment['uri'])

    # Keep the sequences in memory so we can look them up when processing the proteins
    protein_sequences = {}
    # Keep track of unique samples
    samples = {}

    '''
    for resource in dp.resources:
        data_type = resource.descriptor.get('odo-dt:dataType', None)
        if data_type is not None:
            data_type_id = data_type.get('@id', None)
            if data_type_id == 'http://ocean-data.org/schema/data-type/v1.0/ProteinSpectralCounts':
                protein_data_url = mc.presigned_get_object('ingested', 'Metzyme/' + resource.descriptor['path'], expires=timedelta(minutes=1))
                print('proteins data found', protein_data_url)
                table = Table(protein_data_url, schema=resource.schema)
                row_count = 0
                try:
                    for keyed_row in table.iter(keyed=True):
                        row_count += 1
                        for field_name, field_value in keyed_row.items():
                            field = table.schema.get_field(field_name)
                            if (None is field or 'rdfType' not in field.descriptor or field.descriptor['rdfType'] not in ontology_mappings[data_type_id]):
                                # We don't care about fields we don't know anything about
                                continue
                            # Process the field value based on its field definition
                            processed_value = oceanproteinportal.datapackage.processField(value=field_value, descriptor=field.descriptor, field_type=field.descriptor['rdfType'])

                            # What is it?

                except Exception as e:
                    raise e

            elif data_type_id == 'http://ocean-data.org/schema/data-type/v1.0/FASTA-ProteinIdentifications':
                print('proteins IDs found')
                fasta_url = mc.presigned_get_object('ingested', 'Metzyme/' + resource.descriptor['path'], expires=timedelta(minutes=1))
                print(fasta_url)

                #for record in SeqIO.parse(fh, "fasta"):
                #    #protein_sequences[record.id] = record.seq
                #    print(record.id, record.seq)

            elif data_type_id == 'http://ocean-data.org/schema/data-type/v1.0/PeptideSpectralCounts':
                print('peptide data found')
                #table = Table(resource.descriptor['path'], schema=resource.descriptor['schema'])
    '''

if __name__ == '__main__':
    main()

