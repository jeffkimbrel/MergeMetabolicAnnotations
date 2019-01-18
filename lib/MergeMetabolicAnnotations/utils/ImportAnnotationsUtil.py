import os
import datetime
import logging
import json

#from GenomeAnnotationAPI.GenomeAnnotationAPIClient import GenomeAnnotationAPI
from installed_clients.GenomeAnnotationAPIClient import GenomeAnnotationAPI

class ImportAnnotationsUtil:

    workdir = 'tmp/work'
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

        genome_dict = self.genome_api.get_genome_v1({"genomes": [{"ref": genome_ref}], 'downgrade': 0})[
                "genomes"][0]['data']
        #logging.info(genome_dict.keys())

        return genome_dict

    def get_ontology_dict(self, ontology): # getting only the 'ontology' param from the UI

        fullpath = "/kb/module/test/test_data/"
        ontology_path = fullpath + "/" + self.ontology_lookup[ontology]
        ontology_dict_raw = json.loads(open(ontology_path, "r").read() )

        ontology_dict = {}

        for entry in ontology_dict_raw['term_hash']:

            id = ontology_dict_raw['term_hash'][entry]['id']
            name = ontology_dict_raw['term_hash'][entry]['name']
            ontology_dict[id] = name

        return(ontology_dict)

    def get_annotations_file(self, params):
        filename = '/kb/module/test/test_data/PT3_2.Spades.prokka.kegg.txt'  # docker path
        return [line.strip() for line in open(filename)]

    def annotations_to_genes(self, annotations_raw):
        #annotations_raw = [line.strip() for line in open(args.annotations)]
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

    def update_genome(self, genome_dict):
        for feature in genome_dict['features']:

            geneID = feature['id']
            logging.info(geneID)

            #if geneID in self.genes:



                # if genes[geneID].hasValidAnnotations() == True:
                #
                #     # create some things if they don't exist
                #     if 'ontology_terms' not in feature:
                #         feature['ontology_terms'] = {}
                #
                #     if args.namespace not in feature['ontology_terms']:
                #         feature['ontology_terms'][args.namespace] = {}
                #
                #     for annotation in genes[geneID].ontologyChecked:
                #         if annotation['valid'] == 1:
                #
                #             # add to ontologies present
                #             if args.namespace not in genome_dict['ontologies_present']:
                #                 genome_dict['ontologies_present'][args.namespace] = {}
                #
                #             if annotation['id'] not in genome_dict['ontologies_present'][args.namespace]:
                #                 genome_dict['ontologies_present'][args.namespace][annotation['id']] = annotation['name']
                #
                #             if annotation['id'] not in feature['ontology_terms'][args.namespace]:
                #                 feature['ontology_terms'][args.namespace][annotation['id']] = [current_ontology_event]
                #             else:
                #                 feature['ontology_terms'][args.namespace][annotation['id']].append(current_ontology_event)


    def run(self, params):

        # read and prepare objects/files
        self.validate()
        genome_dict = self.get_genome(params['genome'])
        ontology_dict = self.get_ontology_dict(params['ontology'])
        annotations = self.get_annotations_file(params)

        self.annotations_to_genes(annotations)
        self.add_ontology_event(genome_dict, params['ontology'])

        # process
        for gene in self.genes:
            self.genes[gene].validateGeneID(genome_dict)
            self.genes[gene].validateAnnotationID(ontology_dict)

        self.update_genome(genome_dict)




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
