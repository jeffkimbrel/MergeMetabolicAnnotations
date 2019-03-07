import os
import datetime
import logging
import json
import uuid
import yaml
import re

from installed_clients.GenomeAnnotationAPIClient import GenomeAnnotationAPI
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.WorkspaceClient import Workspace as Workspace
from installed_clients.KBaseReportClient import KBaseReport

class CompareAnnotationsUtil:

    workdir     = 'tmp/work/'
    staging_dir = "/staging/"
    datadir     = "/kb/module/data/"

    translation_locations = {'ec'      : 'EBI_EC.ModelSEED.json',
                             'keggro'  : 'KEGG_RXN.ModelSEED.json',
                             'keggko'  : 'KEGG_KO.ModelSEED.json',
                             'metacyc' : 'Metacyc_RXN.ModelSEED.json',
                             'SSO'     : 'SSO.ModelSEED.json'}

    def __init__(self, config):
        os.makedirs(self.workdir, exist_ok = True)
        self.config       = config
        self.timestamp    = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch      = config['scratch']
        self.genome_api   = GenomeAnnotationAPI(self.callback_url)
        self.dfu          = DataFileUtil(self.callback_url)
        self.gfu          = GenomeFileUtil(self.callback_url)
        self.kbr          = KBaseReport(self.callback_url)
        self.ws_client    = Workspace(config["workspace-url"])

        self.events = {}
        self.genes  = {}
        self.rxns   = {}

    def get_genome(self, genome_ref):
        self.genome = self.genome_api.get_genome_v1({"genomes": [{"ref": genome_ref}], 'downgrade': 0})["genomes"][0]['data']

    def get_ontology_events(self, genome_dict):
        if 'ontology_events' in genome_dict:
            for event, ontology in enumerate(genome_dict['ontology_events']):
                self.events[event] = {}
                for term in ontology:
                    self.events[event][term] = ontology[term]
        else:
            logging.info("No ontology events in this genome!")

    def get_translations(self):
        self.translations = {}

        for type in self.translation_locations:
            translations_path = self.datadir + "/" + self.translation_locations[type]
            ontology_translations = json.loads(open(translations_path, "r").read() )
            self.translations[type] = {}

            for term in ontology_translations['translation']:
                for entry in ontology_translations['translation'][term]['equiv_terms']:
                    rxn = entry['equiv_term']
                    self.translations[type][term] = rxn

    def get_genes_and_terms(self):
        for feature in self.genome['features']:
            gene = feature['id']
            if 'ontology_terms' in feature:
                for type in feature['ontology_terms']:
                    for term in feature['ontology_terms'][type]:
                        for ontology_event in feature['ontology_terms'][type][term]:
                            if gene in self.genes:
                                self.genes[gene].add(ontology_event, type, term)
                            else:
                                self.genes[gene] = Gene(gene)
                                self.genes[gene].add(ontology_event, type, term)

    def search_for_ec(self, line):
        ecList = re.findall(r"\(*[0-9]+\.[0-9\-]+\.[0-9\-]+\.[0-9\-]+", line)
        return(ecList)

    def translate_to_rxns(self, getECs = True):
        self.rxns['None'] = RXN('None')

        for gene in self.genes:
            for event in self.genes[gene].annotations:
                term = event['term']
                type = event['type']

                if type == 'SSO':
                    term = self.genome['ontologies_present']['SSO'][term]

                if type == 'metacyc':
                    if term.startswith("META:"):
                        term = term.replace('META:', '')

                if term in self.translations[type]:
                    rxn = self.translations[type][term]

                    if rxn != None:
                        if rxn in self.rxns:
                            self.rxns[rxn].add(event['ontology_event'], gene, type, term)
                        else:
                            self.rxns[rxn] = RXN(rxn)
                            self.rxns[rxn].add(event['ontology_event'], gene, type, term)
                    else:
                        self.rxns['None'].add(event['ontology_event'], gene, type, term)

                else:
                    self.rxns['None'].add(event['ontology_event'], gene, type, term)

                if getECs: #extract ECs from SSO terms
                    if type == "SSO":
                        # ECs
                        ecList = self.search_for_ec(term)
                        if len(ecList) > 0:
                            for ec in ecList:
                                if ec in self.translations['ec']:
                                    rxn = self.translations['ec'][ec]

                                    if rxn != None:
                                        if rxn in self.rxns:
                                            self.rxns[rxn].add(event['ontology_event'], gene, type, ec)
                                        else:
                                            self.rxns[rxn] = RXN(rxn)
                                            self.rxns[rxn].add(event['ontology_event'], gene, type, ec)
                                    else:
                                        self.rxns['None'].add(event['ontology_event'], gene, type, ec)

                                else:
                                    self.rxns['None'].add(event['ontology_event'], gene, type, ec)

    def html_summary(self, params):
        summary = {}

        for gene in self.genes:
            for event in self.genes[gene].annotations:
                term = event['term']
                event = event['ontology_event']

                if event not in summary:
                    summary[event] = {}
                if 'gene' not in summary[event]:
                    summary[event]['gene'] = []
                if 'term' not in summary[event]:
                    summary[event]['term'] = []

                summary[event]['gene'].append(gene)
                summary[event]['term'].append(term)

                summary[event]['gene'] = list(set(summary[event]['gene']))
                summary[event]['term'] = list(set(summary[event]['term']))

        for rxn in self.rxns:
            for event in self.rxns[rxn].translations:
                event = event['ontology_event']
                if event not in summary:
                    summary[event] = {}
                if 'rxn' not in summary[event]:
                    summary[event]['rxn'] = []

                if rxn != "None":
                    summary[event]['rxn'].append(rxn)

                summary[event]['rxn'] = list(set(summary[event]['rxn']))

        output_html_files = list()

        # Make report directory and copy over files
        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        os.mkdir(output_directory)
        result_file_path = os.path.join(output_directory, 'compare_annotations_summary.html')

        # make html
        table_lines = []
        table_lines.append(f'<h2>Compare Annotations</h2>')
        table_lines.append(f'<h3>Summary</h3>')
        table_lines.append('<table cellspacing="0" cellpadding="3" border="1"><tr><th>EVENT</th><th>DESCRIPTION</th><th>TYPE</th><th>GENES</th><th>TERMS</th><th>RXNS</th></tr>')
        for event in sorted(summary.keys()):
            description = self.events[event].get('description', self.events[event]['method']) # RAST/PROKKA don't have descriptions, but they have methods
            type = self.events[event]['id']
            genes_list = summary[event]['gene']
            terms_list = summary[event]['term']
            rxns_list = summary[event]['rxn']
            table_lines.append('<tr><td>' + str(event) + '</td><td>' + description + '</td><td>' + type + '</td><td>' + str(len(set(genes_list))) + '</td><td>' + str(len(terms_list)) + '</td><td>' +  str(len(rxns_list)) + '</td></tr>')
        table_lines.append('</table>')

        # Write to file
        with open(result_file_path, 'w') as result_file:
            for line in table_lines:
                result_file.write(line + "\n")

        output_html_files.append(
            {'path'        : output_directory,
             'name'        : os.path.basename(result_file_path),
             'description' : 'HTML report for compare_metabolic_annotations app'})

        report_params = {
            'message'                   : '',
            'html_links'                : output_html_files,
            'direct_html_link_index'    : 0,
            'workspace_name'            : params['workspace_name'],
            'report_object_name'        : f'compare_annotations_{uuid.uuid4()}'}

        output = self.kbr.create_extended_report(report_params)

        return {'report_name' : output['name'],
                'report_ref'  : output['ref']}


    def run(self, ctx, params):

        # collect some metadata
        self.get_genome(params['genome'])
        self.get_ontology_events(self.genome)
        self.get_translations()

        # collect some more data
        self.get_genes_and_terms()
        self.translate_to_rxns(getECs = True)

        # make reports
        report = self.html_summary(params)

        return report

class Gene:
    def __init__(self, id):
        self.id = id
        self.annotations = []

    def add(self, ontology_event, type, term):
        self.annotations.append({"ontology_event" : ontology_event,
                                           "type" : type,
                                           "term" : term
                               })

class RXN:
    def __init__(self, rxn):
        self.rxn = rxn
        self.translations = []

    def add(self, ontology_event, gene, type, term):
        self.translations.append({"ontology_event" : ontology_event,
                                            "gene" : gene,
                                            "type" : type,
                                            "term" : term
                               })
