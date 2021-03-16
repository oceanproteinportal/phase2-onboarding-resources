from Bio import SeqIO
from rdflib import BNode, Graph, Literal, Namespace, RDF, URIRef, plugin
from rdflib.namespace import XSD
from rdflib.serializer import Serializer
import datapackage
import tableschema.exceptions
from tableschema import Table
from datetime import timedelta
import uuid
from minio import Minio
from minio.error import S3Error

def main():

    mc = Minio(
        "example.com",
        access_key="xxxx",
        secret_key="xxxx",
    )

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

    # Generate datasetId
    datasetId = str(uuid.uuid3(uuid.NAMESPACE_DNS, dp.descriptor['name'] + '_ver.' + dp.descriptor.get('version', 'noversion')))

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
    g = Graph()
    g.bind("geo", geo)
    g.bind("opp", opp)
    g.bind("rdf", RDF)
    g.bind("rdfs", rdfs)
    g.bind("schema", schema)
    g.bind("sosa", sosa)
    g.bind("ssn", ssn)
    g.bind("ssn-system", ssn_system)
    g.bind("wikidata", wikidata)
    g.bind("xsd", XSD)

    # Define URI space for OPP subjects
    id_uri_prefix = "http://id.oceanproteinportal.org/"

    dataset = URIRef(id_uri_prefix + "dataset/" + datasetId)

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

    #Store any cruises
    cruises = {}
    if 'odo:hasDeployment' in dp.descriptor:
        for deployment in dp.descriptor['odo:hasDeployment']:
            cruises[deployment['name']] = URIRef(deployment['uri'])

        print(cruises)

    print(g.serialize(format='ttl').decode("utf-8"))

    # Keep the sequences in memory so we can look them up when processing the proteins
    protein_sequences = {}

    for resource in dp.resources:
        data_type = resource.descriptor.get('odo-dt:dataType', None)
        if data_type is not None:
            data_type_id = data_type.get('@id', None)
            if data_type_id == 'http://ocean-data.org/schema/data-type/v1.0/ProteinSpectralCounts':
                print('proteins data found')
                #table = Table(resource.descriptor['path'], schema=resource.descriptor['schema'])

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


if __name__ == '__main__':
    main()

