import os
import datetime
import logging
import json
from installed_clients.GenomeAnnotationAPIClient import GenomeAnnotationAPI

class ImportAnnotationsUtil:

    workdir = 'tmp/work'
    datadir = "/kb/module/data/"
    ontology_lookup = {
        "ec": "KBaseOntology/ec_ontology",
        "keggko": "KEGG_KO_ontologyDictionary.json",
        "keggro": "KBaseOntology/keggro_ontology",
        "metacyc": "KBaseOntology/metacyc_ontology",
        "modelseed": "KBaseOntology/modelseed_ontology"
    }

    def __init__(self, config):
        os.makedirs(self.workdir, exist_ok = True)
        self.config = config
        self.timestamp = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
        self.genes = {}
        self.callback_url = config['SDK_CALLBACK_URL']
        self.genome_api = GenomeAnnotationAPI(self.callback_url)

    def validate(self):
        pass

    def get_genome(self, genome_ref):
        # TODO: create ontology_event and ontologies_present if they don't exist?

        self.genome_full = self.genome_api.get_genome_v1({"genomes": [{"ref": genome_ref}], 'downgrade': 0})[
                "genomes"][0]
        genome_dict = self.genome_full['data']

        return genome_dict

    def get_ontology_dict(self, ontology): # getting only the 'ontology' param from the UI

        ontology_path = self.datadir + "/" + self.ontology_lookup[ontology]
        ontology_dict_raw = json.loads(open(ontology_path, "r").read() )

        ontology_dict = {}

        for entry in ontology_dict_raw['term_hash']:

            id = ontology_dict_raw['term_hash'][entry]['id']
            name = ontology_dict_raw['term_hash'][entry]['name']
            ontology_dict[id] = name

        return(ontology_dict)

    def get_annotations_file(self, params):
        if params['debug'] == True:
            filename = '/kb/module/test/test_data/marinobacter.prokka.kegg.txt'  # docker path
        else:
            download_staging_file_params = {
                'staging_file_subdir_path': params.get('annotation_file')
            }
            filename = self.dfu.download_staging_file(download_staging_file_params).get('copy_file_path')

        return [line.strip() for line in open(filename)]

    def annotations_to_genes(self, annotations_raw):

        for line in annotations_raw:
            if not line.startswith('#'): # ignore comment lines
                elements = line.split("\t") # can add commas here as well for .csv

                if elements[0] != "": # ignore blank lines
                    geneID = elements[0]
                    annotation = ""

                    # make Gene class for geneID
                    if geneID not in self.genes:
                        self.genes[geneID] = Gene(geneID)

                    # and add the (not yet validated) annotations, columns above 2 are ignored
                    if len(elements) > 1:
                        annotation = elements[1]
                        self.genes[geneID].addAnnotation(annotation)

    def add_ontology_event(self, genome_dict, ontology):
        genome_dict['ontology_events'].append(
            {
                "id"             : ontology,
                "method"         : "TEST",
                "method_version" : "TEST",
                "ontology_ref"   : "TEST",
                "timestamp"      : self.timestamp
            }
        )

    def update_genome(self, genome_dict, ontology):
        for feature in genome_dict['features']:

            geneID = feature['id']

            if geneID in self.genes:

                if self.genes[geneID].hasValidAnnotations() == True:

                    # create some things if they don't exist
                    if 'ontology_terms' not in feature:
                        feature['ontology_terms'] = {}

                    if ontology not in feature['ontology_terms']:
                        feature['ontology_terms'][ontology] = {}

                    for annotation in self.genes[geneID].ontologyChecked:
                        if annotation['valid'] == 1:

                            # add to ontologies present
                            if ontology not in genome_dict['ontologies_present']:
                                genome_dict['ontologies_present'][ontology] = {}

                            if annotation['id'] not in genome_dict['ontologies_present'][ontology]:
                                genome_dict['ontologies_present'][ontology][annotation['id']] = annotation['name']

                            if annotation['id'] not in feature['ontology_terms'][ontology]:
                                feature['ontology_terms'][ontology][annotation['id']] = [self.current_ontology_event]
                            else:
                                feature['ontology_terms'][ontology][annotation['id']].append(self.current_ontology_event)

    def summarize(self):

        validGeneCount = 0
        invalidGenes = []
        validOntologyTermCount = 0
        invalidOntologyTerms = []

        for gene in self.genes:
            if self.genes[gene].valid == 1:
                validGeneCount += 1
            elif self.genes[gene].valid == 0:
                invalidGenes.append(self.genes[gene].id)

            for annotation in self.genes[gene].ontologyChecked:
                if annotation['valid'] == 1:
                    validOntologyTermCount += 1
                elif annotation['valid'] == 0:
                    invalidOntologyTerms.append(annotation['id'])

        logging.info(str(validGeneCount) + "\t" + str(len(invalidGenes)) + "\t" + str(invalidGenes))
        logging.info(str(validOntologyTermCount) + "\t" + str(len(invalidOntologyTerms)) + "\t" + str(invalidOntologyTerms))

    def run(self, ctx, params):

        # read and prepare objects/files
        self.validate()
        genome_dict = self.get_genome(params['genome'])
        ontology_dict = self.get_ontology_dict(params['ontology'])
        annotations = self.get_annotations_file(params)

        self.annotations_to_genes(annotations)
        self.add_ontology_event(genome_dict, params['ontology'])
        self.current_ontology_event = len(genome_dict['ontology_events']) - 1

        # process
        for gene in self.genes:
            self.genes[gene].validateGeneID(genome_dict)
            self.genes[gene].validateAnnotationID(ontology_dict)

        self.update_genome(genome_dict, params['ontology'])

        self.summarize()

        # save genome to new object
        self.genome_full['data'] = genome_dict

        # TODO - ask for these from user, or pull from original genome
        # and what is id?
        for item in ['id', 'scientific_name', 'domain', 'genetic_code']:
            if item not in self.genome_full:
                self.genome_full[item] = "unknown"

                if item == 'genetic_code':
                    self.genome_full[item] = 11

        #logging.info(self.genome_full.keys())
        #logging.info(self.genome_full['data'].keys())

        # with open("/kb/module/work/genome_full.json", 'w') as outfile1:
        #     json.dump(self.genome_full, outfile1, indent = 2)
        #
        # with open("/kb/module/work/genome_dict.json", 'w') as outfile2:
        #     json.dump(genome_dict, outfile2, indent = 2)


        prov = ctx.provenance()
        info = self.genome_api.save_one_genome_v1({'workspace': params['workspace_name'],
                                      'name': params['output_name'],
                                      'data': self.genome_full, 'provenance': prov})['info']

        genome_ref = str(info[6]) + '/' + str(info[0]) + '/' + str(info[4])
        logging.info(genome_ref)

        return {}

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

    def validateAnnotationID(self, ontology_dict):
        for id in self.annotations:
            valid = 0
            if id in ontology_dict:
                valid = 1
                name = ontology_dict[id]

            else:
                name = ""

            # adds valid and not valid annotations for recordkeeping
            ontologyCheck = {"id"    : id,
                             "name"  : name,
                             "valid" : valid
                            }

            self.ontologyChecked.append(ontologyCheck)

    def hasValidAnnotations(self):
        valid = False
        if len(self.ontologyChecked) > 0:
            for annotation in self.ontologyChecked:
                if annotation['valid'] == 1:
                    valid = True
        return(valid)
