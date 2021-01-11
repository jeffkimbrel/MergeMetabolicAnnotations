import yaml
import json
import logging
import pandas as pd
import os
import re


class Gene:
    '''
    This holds uploaded annotation data, where each gene/locus_tag is an object
    '''

    def __init__(self, id):
        self.id = id
        self.valid = 0
        self.annotations = []
        self.ontology_checked = []
        self.feature = None

    def add_annotation(self, annotation):
        '''
        Add annotation term to list, and make unique
        '''
        self.annotations.append(annotation)
        self.annotations = list(set(self.annotations))

    def validate_gene_ID(self, genome):
        '''
        verify that the locus_tag is found in the genome.
        '''
        for feature in genome['features']:
            if feature['id'] == self.id:
                self.valid = 1
                self.feature = feature['id']
            else:  # this changes the import gene id to match the feature id if it is found as an alias
                for alias in feature['aliases']:
                    if self.id == alias[1]:
                        self.feature = feature['id']
                        self.valid = 1

    def validate_annotation_ID(self, ontology_dict, ontology):
        for id in self.annotations:
            name = ""
            valid = 0
            id_final = id

            if id in ontology_dict:
                valid = 1
                name = ontology_dict[id]

            elif ontology == 'META':
                metacyc_id = "META:" + id
                if metacyc_id in ontology_dict:
                    valid = 1
                    name = ontology_dict[metacyc_id]
                    id_final = metacyc_id

            elif ontology == 'GO':
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

            self.ontology_checked.append(ontologyCheck)

    def has_valid_annotations(self):
        valid = False
        if len(self.ontology_checked) > 0:
            for annotation in self.ontology_checked:
                if annotation['valid'] == 1:
                    valid = True
        return(valid)


def validate():
    pass


def get_app_version():
    with open("/kb/module/kbase.yml", 'r') as stream:
        data_loaded = yaml.load(stream)
    return(data_loaded['module-version'])


def get_genome(genome_ref, genome_api):
    '''
    returns just the 'data' portion of the genome object
    '''

    genome = genome_api.get_genome_v1({"genomes": [{"output[“ftrs_not_found”]": genome_ref}],
                                       'downgrade': 0})["genomes"][0]

    return genome['data']


def get_ontology_dict(ontology_type, datadir, ontology_lookup):
    '''
    returns the ontology dictionary selected in the app params. This converts
    everything to upper case to make case insensitive
    key: value = id: name
    '''

    # make case insensitive by converting to upper case
    ontology_type = ontology_type.upper()
    ontology_path = os.path.join(datadir, ontology_lookup[ontology_type])
    ontology_dict_raw = json.loads(open(ontology_path, "r").read())['term_hash']

    ontology_dict = {}

    for entry in ontology_dict_raw:
        id = ontology_dict_raw[entry]['id']

        # Add ontology type if not present

        id = standardize_annotation(id, ontology_type)

        name = ontology_dict_raw[entry]['name']
        ontology_dict[id] = name

    return ontology_dict


def get_annotations_file(params, staging_dir, pass_df=None):
    '''
    Returns the uploaded annotation file as a pandas dictionary with gene names
    in the first column, and ontology terms in the second.  Also adds the
    "ontology:" prefix, if missing

    pass_df option is new accepting a pandas dataframe rather than a file from staging
    '''

    if isinstance(pass_df, pd.DataFrame):
        annotations = pass_df
    else:
        if 'debug' in params and params['debug'] is True:
            annotations_file_path = os.path.join(
                '/kb/module/test/test_data', params.get('annotation_file'))
        else:
            annotations_file_path = os.path.join(staging_dir, params.get('annotation_file'))

        annotations = pd.read_csv(annotations_file_path,
                                  sep='\t',
                                  header=None,
                                  names=['gene', 'term']
                                  )

    # add prefix
    ontology_type = params['ontology'].upper()
    for index, row in annotations.iterrows():
        if pd.isnull(row['term']):
            row['term'] = "NA"

        row['term'] = standardize_annotation(row['term'], params['ontology'].upper())

    return annotations


def annotations_to_genes(annotations, genes):
    '''
    Adds the annotation terms to the gene classes.
    '''
    for index, row in annotations.iterrows():

        # make Gene class for gene
        if row['gene'] not in genes:
            genes[row['gene']] = Gene(row['gene'])

        # and add the (not yet validated) annotations, columns above 2 are ignored
        genes[row['gene']].add_annotation(row['term'])

    return(genes)


def add_ontology_event(genome, params, timestamp, method):
    '''
    Adds the ontology event, creating ontology_events if necessary
    '''
    if 'ontology_events' not in genome:
        genome['ontology_events'] = []

    if 'ontologies_present' not in genome:
        genome['ontologies_present'] = {}

    genome['ontology_events'].append(
        {
            "id": params['ontology'].upper(),
            "method": method,
            "method_version": get_app_version(),
            "description": params['description'],
            # "ontology_ref": "605/6/1",  # this is a hack for now... it is the SSO reference, but this field appears necessary
            "timestamp": timestamp
        }
    )

    return(genome)

##


def update_genome(genome, ontology, genes, current_ontology_event):
    # logging.info(genes.keys())
    for feature in genome['features']:

        featureID = feature['id']

        for gene in genes:
            if genes[gene].has_valid_annotations():
                if featureID == genes[gene].feature:
                    # create some things if they don't exist
                    if 'ontology_terms' not in feature:
                        feature['ontology_terms'] = {}

                    if ontology not in feature['ontology_terms']:
                        feature['ontology_terms'][ontology] = {}

                    for annotation in genes[gene].ontology_checked:
                        if annotation['valid'] == 1:

                            # add to ontologies present
                            if ontology not in genome['ontologies_present']:
                                genome['ontologies_present'][ontology] = {}

                            if annotation['id'] not in genome['ontologies_present'][ontology]:
                                genome['ontologies_present'][ontology][annotation['id']
                                                                       ] = annotation['name']

                            if annotation['id'] not in feature['ontology_terms'][ontology]:
                                feature['ontology_terms'][ontology][annotation['id']] = [
                                    current_ontology_event]
                            else:
                                feature['ontology_terms'][ontology][annotation['id']].append(
                                    current_ontology_event)

    return(genome)


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

        for annotation in genes[gene].ontology_checked:
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


def get_bulk_annotations_file(params, staging_dir):
    '''
    Also adds the "ontology:" prefix, if missing
    '''

    if 'debug' in params and params['debug'] is True:
        annotations_file_path = '/kb/module/test/test_data/' + params.get('annotation_file')
    else:
        annotations_file_path = staging_dir + "/" + params.get('annotation_file')

    df = pd.read_csv(annotations_file_path,
                     sep='\t',
                     header=None,
                     names=['description', 'ontology', 'gene', 'term']
                     )

    # add prefix
    for index, row in df.iterrows():

        if pd.isnull(row['term']):
            row['term'] = "NA"

        row['term'] = standardize_annotation(row['term'], row['ontology'])

    return df


def validate_bulk(bulk_annotations):
    # TODO - had validation check that the first three columns have no null values
    # pd.isnull(df).sum() > 0
    pass


def get_description_ontology_pairs(bulk_annotations):
    return bulk_annotations.groupby(['description', 'ontology']).size().reset_index().rename(columns={0: 'count'})


translation_locations = {'EC': 'EBI_EC.ModelSEED.json',
                         'RO': 'KEGG_RXN.ModelSEED.json',
                         'KO': 'KEGG_KO.ModelSEED.json',
                         'META': 'Metacyc_RXN.ModelSEED.json',
                         'SSO': 'SSO.ModelSEED.json',
                         'GO': 'GO.ModelSEED.json'}

ontology_lookup = {
    "EC": "EBI_EC_ontologyDictionary.json",
    "KO": "KEGG_KO_ontologyDictionary.json",
    "RO": "KEGG_RXN_ontologyDictionary.json",
    "META": "MetaCyc_RXN_ontologyDictionary.json",
    "MSRXN": "ModelSEED_RXN_ontologyDictionary.json",
    "GO": "GO_ontologyDictionary.json",
    "SSO": "SSO_ontologyDictionary.json"
}


def search_for_ec(line):
    ecList = re.findall(r"\(*[0-9]+\.[0-9\-]+\.[0-9\-]+\.[0-9\-]+", line)
    return(ecList)


def get_translations(datadir):
    translations = {}

    for ontology_type in translation_locations:
        translations_path = datadir + "/" + translation_locations[ontology_type]
        ontology_translations = json.loads(open(translations_path, "r").read())
        translations[ontology_type] = {}

        for term in ontology_translations['translation']:

            for entry in ontology_translations['translation'][term]['equiv_terms']:
                if entry['equiv_term'] != None:
                    rxn = entry['equiv_term']
                    # logging.info(rxn)
                    if not rxn.upper().startswith('MSRXN:'):
                        rxn = 'MSRXN:' + rxn
                    # logging.info(rxn)

                    # Add ontology type if not present
                    term = standardize_annotation(term, ontology_type)

                    translations[ontology_type][term] = rxn

    return translations


legacy_codes = {'MODELSEED': 'MSRXN',
                'KEGGKO': 'KO',
                'KEGGRO': 'RO',
                'METACYC': 'META'
                }


def legacy_fix(old):
    '''
    update ontology type with updated terms
    '''
    new = old.upper()
    if new in legacy_codes:
        new = legacy_codes[new]

    return new


def standardize_annotation(annotation, ontology_type):
    '''
    fix annotation to conform to our standards
    '''

    if not annotation.startswith(ontology_type + ':'):
        if annotation.startswith('ec:'):
            annotation = annotation.upper()
        else:
            annotation = ontology_type.upper() + ":" + annotation

    return annotation
