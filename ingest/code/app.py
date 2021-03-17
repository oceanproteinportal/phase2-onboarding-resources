from Bio import SeqIO
from rdflib import BNode, ConjunctiveGraph, Graph, Literal, Namespace, RDF, URIRef, plugin
from rdflib.namespace import XSD
from rdflib.serializer import Serializer
import datapackage
import tableschema.exceptions
from tableschema import Table
from datetime import timedelta
import datetime
import dateutil.parser
import uuid
import unicodedata
import string
import yaml
import pprint
from minio import Minio
from minio.error import S3Error

# GLOBALS
DATE_TIME_FORMAT = "%Y-%m-%dT%H:%M:%S"
VALID_FILENAME_CHARS = "-_.() %s%s" % (string.ascii_letters, string.digits)
FILENAME_CHAR_LIMIT = 255

# Define Namespaces
geo = Namespace("http://www.w3.org/2003/01/geo/wgs84_pos#")
geosparql = Namespace("http://www.opengis.net/ont/geosparql#")
idorg = Namespace("https://registry.identifiers.org/registry/")
opp = Namespace("http://schema.oceanproteinportal.org/v2/")
rdfs = Namespace("http://www.w3.org/2000/01/rdf-schema#")
schema = Namespace("http://schema.org/")
sf = Namespace("http://www.opengis.net/ont/sf#")
sosa = Namespace("http://www.w3.org/ns/sosa/")
ssn = Namespace("http://www.w3.org/ns/ssn/")
ssn_system = Namespace("http://www.w3.org/ns/ssn/systems/")
wikidata = Namespace("http://www.wikidata.org/entity/")


def clean_filename(filename, whitelist=VALID_FILENAME_CHARS, replace=' '):
    """Clean up a potential filename"""

    # replace spaces
    for r in replace:
        filename = filename.replace(r,'_')

    # keep only valid ascii chars
    cleaned_filename = unicodedata.normalize('NFKD', filename).encode('ASCII', 'ignore').decode()

    # keep only whitelisted chars
    cleaned_filename = ''.join(c for c in cleaned_filename if c in whitelist)
    if len(cleaned_filename) > FILENAME_CHAR_LIMIT:
        print("Warning, filename truncated because it was over {}. Filenames may no longer be unique".format(FILENAME_CHAR_LIMIT))
    return cleaned_filename[:FILENAME_CHAR_LIMIT]

def getOntologyMappings(config_file):
    """Get the mappings from the data templates to the ontology"""

    with open(config_file, 'r') as yamlfile:
        mappings = yaml.load(yamlfile, Loader=yaml.FullLoader)

    # Allow the software to lookup fields by the default column name or the ontology definition
    for data_type in list(mappings):
        data_type_def = mappings[data_type]
        for default_column_name in list(data_type_def):
            mappings[data_type][data_type_def[default_column_name]['class']] = data_type_def[default_column_name]
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

def addOppIdentifierToGraph(graph, uri, rdf_type, scheme, value):
    """Consistently construct OPP Identifiers"""
    identifier = BNode(uri)
    graph.add((identifier, RDF.type, rdf_type))
    graph.add((identifier, opp.identifierScheme, scheme))
    graph.add((identifier, opp.identifierValue, Literal(value, datatype=XSD.token)))
    return identifier

def generateProteinResource(datasetID, sampleID, proteinID):
    """Consistently generate protein resource URIs because they need to be calculated from multiple places in the data"""
    return BNode("Protein-{protein}_Sample-{sample}_Dataset={dataset}".format(dataset=datasetID, sample=sampleID, protein=proteinID))

def writeDataFile(graph, rdf_format, directory, filename_prefix='data', row_start='0', row_end='end'):
    dest = "{dir}/{ftype}_{start}_{end}.nq".format(dir=directory, ftype=filename_prefix, start=row_start, end=row_end)
    graph.serialize(destination=dest, format=rdf_format)
    print('Wrote the ', filename_prefix, ' file: ', dest)

def main():
    """Dome some work"""

    # TODO - remove the PrettyPrinter
    pp = pprint.PrettyPrinter(indent=2)

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
    dp_url = mc.presigned_get_object('ingested', CONFIG['dataset']['datapackage'], expires=timedelta(minutes=1))
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

    # Add namespaces to graph
    cg.bind("geo", geo)
    cg.bind("ggeosparql", geosparql)
    cg.bind("idorg", idorg)
    cg.bind("opp", opp)
    cg.bind("rdf", RDF)
    cg.bind("rdfs", rdfs)
    cg.bind("schema", schema)
    cg.bind("sf", sf)
    cg.bind("sosa", sosa)
    cg.bind("ssn", ssn)
    cg.bind("ssn-system", ssn_system)
    cg.bind("wikidata", wikidata)
    cg.bind("xsd", XSD)

    # Define URI space for OPP subjects
    id_uri_prefix = "http://id.oceanproteinportal.org/"
    odo_dt = "http://ocean-data.org/schema/data-type/"

    # Generate dataset_id & define a Dataset resource
    dataset_id = str(uuid.uuid3(uuid.NAMESPACE_DNS, dp.descriptor['name'] + '_ver.' + dp.descriptor.get('version', 'noversion')))
    dataset = URIRef(id_uri_prefix + "dataset/" + dataset_id)

    # Define a Named Graph
    g = Graph(cg.store, dataset)

    # Construct the Dataset
    g.add((dataset, RDF.type, schema.Dataset))
    # Dataset Identifier
    datasetIdScheme = BNode("datasetIdScheme-{s}".format(s=dataset_id))
    g.add((datasetIdScheme, RDF.type, opp.DatasetIdentfierScheme))
    datasetIdentifier = addOppIdentifierToGraph(graph=g, uri=dataset_id,
        rdf_type=opp.Identifier, scheme=datasetIdScheme, value=dataset_id)
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
            orcid = addOppIdentifierToGraph(graph=g, uri="orcid-{s}".format(s=contributor['orcid']),
                rdf_type=opp.Identifier, scheme=wikidata.Q51044, value=contributor['orcid'])
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

    # TODO - At this point we really don't need these triples in memory anymore.
    # 1. Can we store them somewhere to be used later and then clear the graph?
    # Store the dataset RDF
    cg.serialize(destination='./rdf/dataset.nq', format='nquads')
    # Clear the in-memory graph
    g.remove((None, None, None))
    print('Wrote the dataset file')

    #Store any cruises
    cruises = {}
    if 'odo:hasDeployment' in dp.descriptor:
        for deployment in dp.descriptor['odo:hasDeployment']:
            cruises[deployment['name']] = URIRef(deployment['uri'])

    # Keep the sequences in memory so we can look them up when processing the proteins
    protein_sequences = {}
    # Keep track of unique samples
    samples = {}
    filter_sizes = {}
    proteins = {}
    chunk_rows = CONFIG['dataset']['chunk_rows']

    datapackage_filepath_prefix = CONFIG['dataset'].get('path_prefix', '')
    for resource in dp.resources:
        data_type = resource.descriptor.get('odo-dt:dataType', None)
        if data_type is not None:
            data_type_id = data_type.get('@id', None)
            if data_type_id == 'http://ocean-data.org/schema/data-type/v1.0/ProteinSpectralCounts':
                protein_data_url = mc.presigned_get_object('ingested', datapackage_filepath_prefix + resource.descriptor['path'], expires=timedelta(minutes=1))
                #print('proteins data found', protein_data_url)
                table = Table(protein_data_url, schema=resource.schema)
                row_count = 0

                # Look at the protein data file schema
                #pp.pprint(resource.schema.descriptor)
                #print()

                try:
                    for keyed_row in table.iter(keyed=True):
                        row_count += 1

                        # Is it time to clear memory and write to disk?
                        if row_count % chunk_rows == 0:
                            writeDataFile(graph=cg, rdf_format='nquads', directory='./rdf', filename_prefix='proteins', row_start=str(row_count - chunk_rows), row_end=str(row_count))
                            # Clear the in-memory graph
                            g.remove((None, None, None))

                        # Start with a fresh row
                        row = {}

                        # Process each known field
                        for field_name, field_value in keyed_row.items():
                            field = table.schema.get_field(field_name)
                            if (None is field or 'rdfType' not in field.descriptor or field.descriptor['rdfType'] not in ontology_mappings[data_type_id]):
                                # We don't care about fields we don't know anything about
                                continue
                            else:
                                # Process the field value based on its field definition
                                row[field.descriptor['rdfType']] = processField(value=field_value, descriptor=field.descriptor, field_type=field.descriptor['rdfType'])

                        # Look at the protein data row
                        #pp.pprint(keyed_row)
                        #print()
                        # Look at the 'ontologized' protein data row
                        #pp.pprint(row)
                        #print()
                        #break

                        # Have we seen this Sample before?
                        sample = None
                        sample_id = row['http://ocean-data.org/schema/data-type/v1.0/SampleIdentifier']
                        if sample_id not in samples:
                            # Sample
                            sample = BNode("Sample-{s}".format(s=sample_id))
                            # Add this sample to our lookup dict
                            samples[sample_id] = sample
                            g.add((dataset, opp.storesResultsForSample, sample))
                            g.add((sample, RDF.type, sosa.Sample))
                            sampleIdentifier = addOppIdentifierToGraph(graph=g, uri="Sample-{s}-Identifier".format(s=sample_id),
                                rdf_type=opp.SampleIdentifier, scheme=datasetIdScheme, value=sample_id)
                            g.add((sample, opp.identifier, sampleIdentifier))
                            # Feature of Interest
                            featureOfInterest = BNode("Sample-{s}-FeatureOfInterest".format(s=sample_id))
                            g.add((sample, sosa.isSampleOf, featureOfInterest))
                            g.add((featureOfInterest, RDF.type, sosa.FeatureOfInterest))
                            lat = row['http://ocean-data.org/schema/data-type/v1.0/LatitudeDecimalDegrees']
                            lon = row['http://ocean-data.org/schema/data-type/v1.0/LongitudeDecimalDegrees']
                            depth = row['http://ocean-data.org/schema/data-type/v1.0/DepthMeters']
                            featureName = "Lon: {lon:f}, Lat: {lat:f}".format(lon=lon, lat=lat)
                            g.add((featureOfInterest, geo.lat, Literal(lat, datatype=XSD.float)))
                            g.add((featureOfInterest, geo.lon, Literal(lon, datatype=XSD.float)))
                            # Setup the geometry for spatial query
                            geometry = BNode("Sample-{s}-FeatureOfInterest-Geometry".format(s=sample_id))
                            g.add((featureOfInterest, geo.hasGeometry, geometry))
                            g.add((geometry, RDF.type, sf.Point))
                            g.add((geometry, geosparql.asWKT, Literal("POINT({lon:f} {lat:f})".format(lon=lon, lat=lat), datatype=geosparql.WKTLiteral)))
                            g.add((featureOfInterest, opp.depth, Literal(depth, datatype=XSD.float)))
                            # Station
                            station = BNode("Sample-{s}-Station".format(s=sample_id))
                            g.add((featureOfInterest, opp.inVicinityOfStation, station))
                            g.add((station, RDF.type, opp.Station))
                            stationName = row.get('http://ocean-data.org/schema/data-type/v1.0/CruiseStationIdentifier', featureName)
                            g.add((station, opp.stationName, Literal(stationName, datatype=XSD.string)))
                            # Sampling
                            sampling = BNode("Sample-{s}-Sampling".format(s=sample_id))
                            g.add((sample, sosa.isResultOf, sampling))
                            g.add((sampling, RDF.type, sosa.Sampling))
                            g.add((sampling, sosa.hasFeatureOfInterest, featureOfInterest))
                            # Result Time
                            if 'http://ocean-data.org/schema/data-type/v1.0/ISODateTimeUTC' in row:
                                resultTime = dateutil.parser.parse(row['http://ocean-data.org/schema/data-type/v1.0/ISODateTimeUTC'])
                                resultTime = resultTime.strftime(DATE_TIME_FORMAT)
                                g.add((sampling, sosa.resultTime, Literal(resultTime, datatype=XSD.dateTime)))
                            elif 'http://ocean-data.org/schema/data-type/v1.0/Date' in row:
                                time = row.get('http://ocean-data.org/schema/data-type/v1.0/Time', '00:00:00')
                                resultTime = dateutil.parser.parse(row['http://ocean-data.org/schema/data-type/v1.0/Date'] + 'T' + time)
                                resultTime = resultTime.strftime(DATE_TIME_FORMAT)
                                g.add((sampling, sosa.resultTime, Literal(resultTime, datatype=XSD.dateTime)))
                            # Sampler
                            sampler = BNode("Sample-{s}-Sampler".format(s=sample_id))
                            g.add((sampling, sosa.madeBySampler, sampler))
                            g.add((sampler, RDF.type, sosa.Sampler))
                            # Cruise
                            if 'http://ocean-data.org/schema/data-type/v1.0/CruiseIdentifier' in row:
                                cruise_id = row['http://ocean-data.org/schema/data-type/v1.0/CruiseIdentifier']
                                # Do we know this cruise from the datapackage.json metadata?
                                cruise = None
                                if cruise_id in cruises:
                                    cruise = cruises[cruise_id]
                                else:
                                    cruise = BNode("Cruise-{s}".format(s=cruise_id))
                                    #Make sure we reuse our cruise(URIRef) within this dataset
                                    cruises[cruise_id] = cruise
                                g.add((sampler, sosa.hasDeployment, cruise))
                                g.add((cruise, RDF.type, opp.Cruise))
                                g.add((cruise, opp.cruiseID, Literal(cruise_id, datatype=XSD.string)))
                            # FilterSize
                            min_fs = row.get('http://ocean-data.org/schema/data-type/v1.0/MinFilterSizeInMicrons', 'unspecified minimum')
                            max_fs = row.get('http://ocean-data.org/schema/data-type/v1.0/MaxFilterSizeInMicrons', 'unspecified maximum')
                            filter_size_label = "{min} - {max} microns".format(min=min_fs, max=max_fs)
                            # Have we seen this filter size before?
                            if filter_size_label in filter_sizes:
                                filter_size = filter_sizes[filter_size_label]
                            else:
                                filter_size = BNode("FilterSize-{s}".format(s=clean_filename(filter_size_label)))
                                #Make sure we reuse our filter-size(URIRef) within this dataset
                                filter_sizes[filter_size_label] = filter_size
                            g.add((sampler, opp.hasFilterSize, filter_size))
                            g.add((filter_size, RDF.type, opp.FilterSize))
                            g.add((filter_size, opp.minimumFilterSize, Literal(min_fs, datatype=XSD.float)))
                            g.add((filter_size, opp.maximumFilterSize, Literal(max_fs, datatype=XSD.float)))
                            g.add((sample, ssn.hasProperty, filter_size))
                            # The data doesn't describe a single filter size. This is a placeholder for the future.
                            # Make sure you add this to the 'filter_size_label' above
                            #g.add((filter_size, schema:value, Literal(, datatype=XSD.float)))


                            # The data doesn't describe procedures. This is a placeholder for the future
                            #g.add((sampling, sosa.usedProcedure, ))

                            # serialize the sample
                            sample_file = "./rdf/sample_{id}.nq".format(id=clean_filename(sample_id))
                            cg.serialize(destination=sample_file, format='nquads')
                            # Clear the in-memory graph
                            g.remove((None, None, None))
                            print('Wrote the sample file: ', sample_file)
                        else:
                            # For reusing the sample resource URI if needed in other places
                            sample = samples[sample_id]

                        # TODO - Set ORF_id to 'http://ocean-data.org/schema/data-type/v1.0/ProteinIdentifier' in Metzyme datapackage.json
                        proteinIdentification = None
                        protein_id = row['http://ocean-data.org/schema/data-type/v1.0/ProteinIdentifier']
                        if protein_id not in proteins:
                            proteinIdentification = generateProteinResource(datasetID=dataset_id, sampleID=sample_id, proteinID=protein_id)
                            # Add this protein identification to our lookup dict
                            proteins[protein_id]= proteinIdentification

                            g.add((proteinIdentification, RDF.type, opp.ProteinIdentification))
                            proteinIdentifier = addOppIdentifierToGraph(graph=g, uri="Protein-{s}-Identifier".format(s=protein_id),
                                rdf_type=opp.ProteinIdentifier, scheme=datasetIdScheme, value=protein_id)
                            g.add((proteinIdentification, opp.identifier, proteinIdentifier))
                            # NOTE the 'Protein' class will be instantiated during the FASTA file read so that the sequence can be attached to it
                            # The FASTA processing will lookup the protein identification in the 'proteins' dict to set 'opp:describesProtein'

                            # protein name
                            if 'http://ocean-data.org/schema/data-type/v1.0/IdentifiedProductName' in row:
                                g.add((proteinIdentification, opp.proteinName, Literal(row['http://ocean-data.org/schema/data-type/v1.0/IdentifiedProductName'], datatype=XSD.string)))
                            # molecular weight
                            if 'http://ocean-data.org/schema/data-type/v1.0/MolecularWeightInDaltons' in row:
                                g.add((proteinIdentification, opp.molecularWeight, Literal(row['http://ocean-data.org/schema/data-type/v1.0/MolecularWeightInDaltons'], datatype=XSD.float)))

                            # Process all the associated Identifications - I think this is the correct place.
                            # TODO - confirm that the identification details do not change per row where the same protein is being identified
                            # Enzyme Commission
                            if 'http://ocean-data.org/schema/data-type/v1.0/EnzymeCommissionIdentifier' in row:
                                addOppIdentifierToGraph(graph=g, uri="Protein-{s}-EnzymeCommissionIdentifier".format(s=protein_id),
                                    rdf_type=opp.ProteinIdentifier, scheme=wikidata.Q741108, value=row['http://ocean-data.org/schema/data-type/v1.0/EnzymeCommissionIdentifier'])
                            # Uniprot
                            if 'http://ocean-data.org/schema/data-type/v1.0/UniprotIdentifier' in row:
                                addOppIdentifierToGraph(graph=g, uri="Protein-{s}-UniprotIdentifier".format(s=protein_id),
                                    rdf_type=opp.ProteinIdentifier, scheme=idorg.uniprot, value=row['http://ocean-data.org/schema/data-type/v1.0/UniprotIdentifier'])

                            #TODO - Reasoning to generate the inverse opp:identifierMetadata triple
                            # PFams
                            if 'http://ocean-data.org/schema/data-type/v1.0/PFamsIdentifier' in row:
                                pfams = addOppIdentifierToGraph(graph=g, uri="Protein-{s}-PFamsIdentifier".format(s=protein_id),
                                    rdf_type=opp.ProteinIdentifier, scheme=idorg.pfam, value=row['http://ocean-data.org/schema/data-type/v1.0/PFamsIdentifier'])
                                if 'http://ocean-data.org/schema/data-type/v1.0/PFamsDescription' in row:
                                    pfamsDescription = BNode("PFam_{s}".format(s=protein_id))
                                    g.add((pfamsDescription, opp.annotatesIdentifier, pfams))
                                    g.add((pfamsDescription, RDF.type, opp.PFamsDescription))
                                    g.add((pfamsDescription, schema.name, Literal(row['http://ocean-data.org/schema/data-type/v1.0/PFamsIdentifier'], datatype=XSD.string)))
                                    g.add((pfamsDescription, schema.description, Literal(row['http://ocean-data.org/schema/data-type/v1.0/PFamsDescription'], datatype=XSD.string)))
                            # Kegg
                            if 'http://ocean-data.org/schema/data-type/v1.0/KeggIdentifier' in row:
                                kegg = addOppIdentifierToGraph(graph=g, uri="Protein-{s}-KeggIdentifier".format(s=protein_id),
                                    rdf_type=opp.ProteinIdentifier, scheme=idorg.kegg, value=row['http://ocean-data.org/schema/data-type/v1.0/KeggIdentifier'])
                                kegg_metadata = BNode("Kegg_{s}".format(s=protein_id))
                                has_kegg_metadata = False
                                if 'http://ocean-data.org/schema/data-type/v1.0/KeggName' in row:
                                    has_kegg_metadata = True
                                    g.add((kegg_metadata, schema.name, Literal(row['http://ocean-data.org/schema/data-type/v1.0/KeggName'], datatype=XSD.string)))
                                if 'http://ocean-data.org/schema/data-type/v1.0/KeggDescription' in row:
                                    has_kegg_metadata = True
                                    g.add((kegg_metadata, schema.description, Literal(row['http://ocean-data.org/schema/data-type/v1.0/KeggDescription'], datatype=XSD.string)))
                                if has_kegg_metadata:
                                    g.add((kegg_metadata, opp.annotatesIdentifier, kegg))
                                    g.add((kegg_metadata, RDF.type, opp.KeggPathway))

                            # NCBI
                            if 'http://ocean-data.org/schema/data-type/v1.0/BestNCBITaxonIdentifier' in row:
                                ncbi = addOppIdentifierToGraph(graph=g, uri="Protein-{s}-NCBITaxonIdentifier".format(s=protein_id),
                                    rdf_type=opp.ProteinIdentifier, scheme=idorg.ncbiprotein, value=row['http://ocean-data.org/schema/data-type/v1.0/BestNCBITaxonIdentifier'])
                                if 'http://ocean-data.org/schema/data-type/v1.0/BestNCBITaxonName' in row:
                                    ncbiTaxonName = BNode("NCBITaxonName_{s}".format(s=protein_id))
                                    g.add((ncbiTaxonName, opp.annotatesIdentifier, ncbi))
                                    g.add((ncbiTaxonName, RDF.type, opp.NCBITaxonName))
                                    g.add((ncbiTaxonName, schema.name, Literal(row['http://ocean-data.org/schema/data-type/v1.0/BestNCBITaxonName'], datatype=XSD.string)))
                        else:
                            proteinIdentification = proteins[protein_id]

                        # Determine the sample protein properties based on the data type.
                        # Actually, we know we are ProteinSpectralCount because of the data type check from processing the datapackage resource property 'odo-dt:dataType'
                        proteinTotalSpectralCount = BNode("Protein-{protein}-Sample-{sample}-TotalSpectralCount".format(protein=protein_id, sample=sample_id))
                        g.add((proteinIdentification, opp.sampleProteinProperty, proteinTotalSpectralCount))
                        g.add((proteinTotalSpectralCount, RDF.type, opp.ProteinTotalSpectralCount))
                        g.add((proteinTotalSpectralCount, opp.observationInSample, sample))
                        g.add((proteinTotalSpectralCount, schema.value, Literal(row['http://ocean-data.org/schema/data-type/v1.0/SpectralCount'], datatype=XSD.float)))

                    # Complete the proteins RDF
                    writeDataFile(graph=cg, rdf_format='nquads', directory='./rdf', filename_prefix='proteins', row_start=str(row_count), row_end='end')
                    # Clear the in-memory graph
                    g.remove((None, None, None))

                except Exception as e:
                    raise e

            elif data_type_id == 'http://ocean-data.org/schema/data-type/v1.0/FASTA-ProteinIdentifications':
                fasta_url = mc.presigned_get_object('ingested', datapackage_filepath_prefix + resource.descriptor['path'], expires=timedelta(minutes=1))
                #print('proteins sequences found', fasta_url)

                #for record in SeqIO.parse(fh, "fasta"):
                #    #protein_sequences[record.id] = record.seq
                #    print(record.id, record.seq)

            elif data_type_id == 'http://ocean-data.org/schema/data-type/v1.0/PeptideSpectralCounts':
                peptide_url = mc.presigned_get_object('ingested', datapackage_filepath_prefix + resource.descriptor['path'], expires=timedelta(minutes=1))
                #print('peptide data found', peptide_url)
                #table = Table(resource.descriptor['path'], schema=resource.descriptor['schema'])

if __name__ == '__main__':
    main()

