import yaml
import json
import logging
import pandas as pd


class Gene:
    def __init__(self, id):
        self.id = id
        self.valid = 0
        self.annotations = []
        self.ontologyChecked = []

    def addAnnotation(self, annotation):
        self.annotations.append(annotation)
        self.annotations = list(set(self.annotations))

    def validateGeneID(self, genome_dict):
        for feature in genome_dict['features']:
            if feature['id'] == self.id:
                self.valid = 1

    def validateAnnotationID(self, ontology_dict, ontology):
        for id in self.annotations:
            name = ""
            valid = 0
            id_final = id

            if id in ontology_dict:
                valid = 1
                name = ontology_dict[id]

            elif ontology == 'metacyc':
                metacyc_id = "META:" + id
                if metacyc_id in ontology_dict:
                    valid = 1
                    name = ontology_dict[metacyc_id]
                    id_final = metacyc_id

            elif ontology == 'go':
                go_id = id.replace('GO:', '')
                if go_id in ontology_dict:
                    valid = 1
                    name = ontology_dict[go_id]
                    id_final = go_id

                # adds valid and not valid annotations for recordkeeping
            ontologyCheck = {"id": id_final,
                             "name": name,
                             "valid": valid
                             }

            self.ontologyChecked.append(ontologyCheck)

    def hasValidAnnotations(self):
        valid = False
        if len(self.ontologyChecked) > 0:
            for annotation in self.ontologyChecked:
                if annotation['valid'] == 1:
                    valid = True
        return(valid)


def validate():
    pass


def get_sso_data(sso_to_lookup, ws_client):
    sso_to_lookup = "KBaseOntology/seed_subsystem_ontology"
    sso_ret = ws_client.get_objects([{"ref": sso_to_lookup}])[0]
    sso = sso_ret["data"]
    sso_info = sso_ret["info"]
    sso_ref = str(sso_info[6]) + "/" + str(sso_info[0]) + "/" + str(sso_info[4])
    return sso_ref


def get_app_version():
    with open("/kb/module/kbase.yml", 'r') as stream:
        data_loaded = yaml.load(stream)
    return(data_loaded['module-version'])


def get_genome(genome_ref, genome_api):
    genome_full = genome_api.get_genome_v1(
        {"genomes": [{"ref": genome_ref}], 'downgrade': 0})["genomes"][0]
    genome_dict = genome_full['data']

    return (genome_dict, genome_full)


def get_ontology_dict(ontology, datadir, ontology_lookup):  # getting only the 'ontology' param from the UI

    ontology_path = datadir + "/" + ontology_lookup[ontology]
    ontology_dict_raw = json.loads(open(ontology_path, "r").read())

    ontology_dict = {}

    for entry in ontology_dict_raw['term_hash']:
        id = ontology_dict_raw['term_hash'][entry]['id']
        name = ontology_dict_raw['term_hash'][entry]['name']
        ontology_dict[id] = name

    return(ontology_dict)


def get_annotations_file(params, staging_dir):
    if 'debug' in params and params['debug'] is True:
        annotations_file_path = '/kb/module/test/test_data/' + params.get('annotation_file')

    else:
        annotations_file_path = staging_dir + "/" + params.get('annotation_file')

    annotations = pd.read_csv(annotations_file_path,
                              sep='\t',
                              header=None,
                              names=['gene', 'term']
                              )

    return annotations
    # return [line.strip() for line in open(annotations_file_path)]


def annotations_to_genes(annotations, genes):
    for index, row in annotations.iterrows():

        # make Gene class for gene
        if row['gene'] not in genes:
            genes[row['gene']] = Gene(row['gene'])

        # and add the (not yet validated) annotations, columns above 2 are ignored
        if len(row['term']) > 1:
            genes[row['gene']].addAnnotation(row['term'])

    # for line in annotations_raw:
    #     if not line.startswith('#'):  # ignore comment lines
    #         elements = line.split("\t")  # can add commas here as well for .csv
    #
    #         if elements[0] != "":  # ignore blank lines
    #             geneID = elements[0]
    #             annotation = ""
    #
    #             # make Gene class for geneID
    #             if geneID not in genes:
    #                 genes[geneID] = Gene(geneID)
    #
    #             # and add the (not yet validated) annotations, columns above 2 are ignored
    #             if len(elements) > 1:
    #                 annotation = elements[1]
    #                 genes[geneID].addAnnotation(annotation)

    return(genes)


def add_ontology_event(genome_dict, params, sso_ref, timestamp):

    if 'ontology_events' not in genome_dict:
        genome_dict['ontology_events'] = []

    if 'ontologies_present' not in genome_dict:
        genome_dict['ontologies_present'] = {}

    genome_dict['ontology_events'].append(
        {
            "id": params['ontology'],
            "method": "Import Annotations",
            "method_version": get_app_version(),
            "description": params['description'],
            "ontology_ref": sso_ref,
            "timestamp": timestamp
        }
    )

    return(genome_dict)


def update_genome(genome_dict, ontology, genes, current_ontology_event):
    for feature in genome_dict['features']:

        geneID = feature['id']

        if geneID in genes:

            if genes[geneID].hasValidAnnotations() is True:

                # create some things if they don't exist
                if 'ontology_terms' not in feature:
                    feature['ontology_terms'] = {}

                if ontology not in feature['ontology_terms']:
                    feature['ontology_terms'][ontology] = {}

                for annotation in genes[geneID].ontologyChecked:
                    if annotation['valid'] == 1:

                        # add to ontologies present
                        if ontology not in genome_dict['ontologies_present']:
                            genome_dict['ontologies_present'][ontology] = {}

                        if annotation['id'] not in genome_dict['ontologies_present'][ontology]:
                            genome_dict['ontologies_present'][ontology][annotation['id']
                                                                        ] = annotation['name']

                        if annotation['id'] not in feature['ontology_terms'][ontology]:
                            feature['ontology_terms'][ontology][annotation['id']] = [
                                current_ontology_event]
                        else:
                            feature['ontology_terms'][ontology][annotation['id']].append(
                                current_ontology_event)

    return(genome_dict)


def summarize(params, genes):
    """
    Goes through the Gene objects and returns a summary dict with
    information on valid and invalid genes and terms. This is used to
    generate the html report later.
    """
    validGenes = []
    invalidGenes = []
    validOntologyTerms = []
    invalidOntologyTerms = []

    for gene in genes:
        if genes[gene].valid == 1:
            validGenes.append(genes[gene].id)
        elif genes[gene].valid == 0:
            invalidGenes.append(genes[gene].id)

        for annotation in genes[gene].ontologyChecked:
            if annotation['valid'] == 1:
                validOntologyTerms.append(annotation['id'])
            elif annotation['valid'] == 0:
                invalidOntologyTerms.append(annotation['id'])

    return({
        'valid_genes': validGenes,
        'invalid_genes': invalidGenes,
        'valid_terms': validOntologyTerms,
        'invalid_terms': invalidOntologyTerms
    })

# bulk functions


def process_bulk_file(params):
    if 'debug' in params and params['debug'] is True:
        annotations_file_path = '/kb/module/test/test_data/' + params.get('annotation_file')
    else:
        annotations_file_path = staging_dir + "/" + params.get('annotation_file')

    df = pd.read_csv(annotations_file_path,
                     sep='\t',
                     header=None,
                     names=['description', 'ontology', 'gene', 'term']
                     )

    # TODO - had validation check that the first three columns have no null values
    logging.info(pd.isnull(df).sum() > 0)

    # split by descriptions
    annotations = {}
    for index, row in df.iterrows():

        if row['description'] not in annotations:
            annotations[row['description']] = {'ontology': [], 'genes': {}}

        if row['ontology'] not in annotations[row['description']]['ontology']:
            annotations[row['description']]['ontology'].append(row['ontology'])

        if row['gene'] not in annotations[row['description']]['genes']:
            annotations[row['description']]['genes'][row['gene']] = []

        annotations[row['description']]['genes'][row['gene']].append(row['term'])

    # logging.info(annotations)

    # verify each description has a single ontology
