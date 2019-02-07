import os
import datetime
import logging
import json
from installed_clients.GenomeAnnotationAPIClient import GenomeAnnotationAPI
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.WorkspaceClient import Workspace as Workspace

class ImportAnnotationsUtil:

    workdir = 'tmp/work'
    staging_dir = "/staging/"
    datadir = "/kb/module/data/"

    ontology_lookup = {
        "ec"        : "EBI_EC_ontologyDictionary.json.txt",
        "keggko"    : "KEGG_KO_ontologyDictionary.json",
        "keggro"    : "KEGG_RXN_ontologyDictionary.json.txt",
        "metacyc"   : "MetaCyc_RXN_ontologyDictionary.json.txt",
        "modelseed" : "ModelSEED_RXN_ontologyDictionary.json.txt"
    }

    def __init__(self, config):
        os.makedirs(self.workdir, exist_ok = True)
        self.config = config
        self.timestamp = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
        self.genes = {}
        self.callback_url = config['SDK_CALLBACK_URL']
        self.genome_api = GenomeAnnotationAPI(self.callback_url)
        self.dfu = DataFileUtil(self.callback_url)
        self.gfu = GenomeFileUtil(self.callback_url)
        self.ws_client = Workspace(config["workspace-url"])

    def get_sso_data(self,sso_to_lookup):
        sso_to_lookup = "KBaseOntology/seed_subsystem_ontology"
        sso_ret = \
            self.ws_client.get_objects([{"ref": sso_to_lookup}])[0]
        sso = sso_ret["data"]
        sso_info = sso_ret["info"]
        self.sso_ref = str(sso_info[6]) + "/" + str(sso_info[0]) + "/" + str(sso_info[4])
        return sso

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
        if 'debug' in params and params['debug'] == True:
            annotations_file_path = '/kb/module/test/test_data/Cdiff_genes_reactions.tsv'  # docker path

        else:
            annotations_file_path = self.staging_dir + "/" + params.get('annotation_file')

        return [line.strip() for line in open(annotations_file_path)]

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

        if 'ontology_events' not in genome_dict:
            genome_dict['ontology_events'] = []


        genome_dict['ontology_events'].append(
            {
                "id"             : ontology,
                "method"         : "TEST",
                "method_version" : "TEST",
                "ontology_ref"   : self.sso_ref,
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
        self.get_sso_data(params['ontology'])

        genome_dict = self.get_genome(params['genome'])
        ontology_dict = self.get_ontology_dict(params['ontology'])
        annotations = self.get_annotations_file(params)

        self.annotations_to_genes(annotations)
        self.add_ontology_event(genome_dict, params['ontology'])
        self.current_ontology_event = len(genome_dict['ontology_events']) - 1

        # process
        for gene in self.genes:
            self.genes[gene].validateGeneID(genome_dict)
            self.genes[gene].validateAnnotationID(ontology_dict, params['ontology'])

        self.update_genome(genome_dict, params['ontology'])

        self.summarize()

        # overwrite object with new genome
        self.genome_full['data'] = genome_dict

        prov = ctx.provenance()
        logging.info(params)

        info = self.gfu.save_one_genome({'workspace' : params['workspace_name'],
                                      'name': params['output_name'],
                                      'data': self.genome_full['data'],
                                      'provenance' : prov})['info']

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

            # adds valid and not valid annotations for recordkeeping
            ontologyCheck = {"id"    : id_final,
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
