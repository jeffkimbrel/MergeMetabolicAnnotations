import os
import datetime
import logging
import json
import uuid
import re

from installed_clients.GenomeAnnotationAPIClient import GenomeAnnotationAPI
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.WorkspaceClient import Workspace as Workspace
from installed_clients.KBaseReportClient import KBaseReport


class CompareAnnotationsUtil:

    workdir = 'tmp/work/'
    staging_dir = "/staging/"
    datadir = "/kb/module/data/"

    translation_locations = {'ec': 'EBI_EC.ModelSEED.json',
                             'keggro': 'KEGG_RXN.ModelSEED.json',
                             'keggko': 'KEGG_KO.ModelSEED.json',
                             'metacyc': 'Metacyc_RXN.ModelSEED.json',
                             'SSO': 'SSO.ModelSEED.json',
                             'go': 'GO.ModelSEED.json'}

    def __init__(self, config):
        os.makedirs(self.workdir, exist_ok=True)
        self.config = config
        self.timestamp = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.genome_api = GenomeAnnotationAPI(self.callback_url)
        self.dfu = DataFileUtil(self.callback_url)
        self.gfu = GenomeFileUtil(self.callback_url)
        self.kbr = KBaseReport(self.callback_url)
        self.ws_client = Workspace(config["workspace-url"])

        self.events = {}

    def get_genome(self, genome_ref):
        self.genome = self.genome_api.get_genome_v1({"genomes": [{"ref": genome_ref}], 'downgrade': 0})[
            "genomes"][0]['data']

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
            ontology_translations = json.loads(open(translations_path, "r").read())
            self.translations[type] = {}

            for term in ontology_translations['translation']:
                for entry in ontology_translations['translation'][term]['equiv_terms']:
                    rxn = entry['equiv_term']
                    self.translations[type][term] = rxn

    def search_for_ec(self, line):
        ecList = re.findall(r"\(*[0-9]+\.[0-9\-]+\.[0-9\-]+\.[0-9\-]+", line)
        return(ecList)

    def summarize_gto(self, gto, translations):
        summary = {"genes": {},
                   "terms": {},
                   "rxns": {},
                   "ontology_events": {},
                   "orphan_terms": {}
                   }

        # add ontology events
        for count, oe in enumerate(gto['ontology_events']):
            summary['ontology_events'][count] = oe

        # add gene id to summary
        for feature in gto['features']:
            gene_id = feature['id']
            summary["genes"][gene_id] = {"terms": {},
                                         "rxns": {}
                                         }

            # get ontology term
            if "ontology_terms" in feature:
                for type in feature['ontology_terms']:
                    term_dict = feature['ontology_terms'][type]

                    for term in term_dict:
                        for oe in term_dict[term]:

                            rxn = "none"

                            # get rxn
                            ontology_type = summary['ontology_events'][oe]['id']

                            # fix metacyc terms
                            if ontology_type == 'metacyc':
                                if term.startswith("META:"):
                                    term = term.replace('META:', '')

                            # fix go terms
                            if ontology_type == 'go':
                                if term.startswith("GO:"):
                                    term = term.replace('GO:', '')

                            # fix SSO terms
                            if ontology_type == 'SSO':

                                if term in gto['ontologies_present']['SSO']:
                                    if gto['ontologies_present']['SSO'][term] != 'Unknown':
                                        term = gto['ontologies_present']['SSO'][term]

                            # convert terms to rxns
                            if term in translations[ontology_type]:
                                rxn = translations[ontology_type][term]
                            else:
                                if oe in summary["orphan_terms"]:
                                    summary["orphan_terms"][oe].append(term)
                                    summary["orphan_terms"][oe] = list(
                                        set(summary["orphan_terms"][oe]))
                                else:
                                    summary["orphan_terms"][oe] = [term]

                            # terms
                            if term in summary["genes"][gene_id]['terms']:
                                summary["genes"][gene_id]['terms'][term].append(oe)
                            else:
                                summary["genes"][gene_id]['terms'][term] = [oe]

                            if term in summary['terms']:
                                summary['terms'][term].append(oe)
                                summary['terms'][term] = list(set(summary['terms'][term]))
                            else:
                                summary['terms'][term] = [oe]

                            # rxns
                            if rxn != "none":
                                if rxn in summary["genes"][gene_id]['rxns']:
                                    summary["genes"][gene_id]['rxns'][rxn].append(oe)
                                else:
                                    summary["genes"][gene_id]['rxns'][rxn] = [oe]

                                if rxn in summary['rxns']:
                                    summary['rxns'][rxn].append(oe)
                                    summary['rxns'][rxn] = list(set(summary['rxns'][rxn]))
                                else:
                                    summary['rxns'][rxn] = [oe]

        with open(os.path.join(self.scratch, "summary_dump.json"), 'w') as outfile:
            json.dump(summary, outfile, indent=2)

        return summary

    def html_summary(self, params, summary):

        # convert gto summary for this report
        html_summary_report = {}

        for oe in summary['ontology_events']:
            html_summary_report[oe] = {"gene": [], "term": [], "rxn": []}

        for gene in summary["genes"]:
            for term in summary["genes"][gene]['terms']:
                for oe in summary["genes"][gene]['terms'][term]:
                    html_summary_report[oe]['gene'].append(gene)
                    html_summary_report[oe]['term'].append(term)

                    html_summary_report[oe]['gene'] = list(set(html_summary_report[oe]['gene']))
                    html_summary_report[oe]['term'] = list(set(html_summary_report[oe]['term']))

            for rxn in summary["genes"][gene]['rxns']:
                for oe in summary["genes"][gene]['rxns'][rxn]:
                    html_summary_report[oe]['rxn'].append(rxn)
                    html_summary_report[oe]['gene'].append(gene)

                    html_summary_report[oe]['rxn'] = list(set(html_summary_report[oe]['rxn']))
                    html_summary_report[oe]['gene'] = list(set(html_summary_report[oe]['gene']))

        output_html_files = list()

        # Make report directory and copy over files
        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        os.mkdir(output_directory)
        result_file_path = os.path.join(output_directory, 'compare_annotations_summary.html')

        # make html
        table_lines = []
        table_lines.append(f'<h2>Compare Annotations</h2>')
        table_lines.append(f'<h3>Summary</h3>')
        table_lines.append(
            '<table cellspacing="0" cellpadding="3" border="1"><tr><th>EVENT</th><th>DESCRIPTION</th><th>TYPE</th><th>GENES</th><th>TERMS</th><th>RXNS</th></tr>')
        for event in sorted(html_summary_report.keys()):
            # RAST/PROKKA don't have descriptions, but they have methods
            description = self.events[event].get('description', self.events[event]['method'])
            type = self.events[event]['id']
            genes_list = html_summary_report[event]['gene']
            terms_list = html_summary_report[event]['term']
            rxns_list = html_summary_report[event]['rxn']
            table_lines.append('<tr><td>' + str(event) + '</td><td>' + description + '</td><td>' + type + '</td><td>' + str(
                len(set(genes_list))) + '</td><td>' + str(len(terms_list)) + '</td><td>' + str(len(rxns_list)) + '</td></tr>')
        table_lines.append('</table>')

        # Write to file
        with open(result_file_path, 'w') as result_file:
            for line in table_lines:
                result_file.write(line + "\n")

        output_html_files.append(
            {'path': output_directory,
             'name': os.path.basename(result_file_path),
             'description': 'HTML report for compare_metabolic_annotations app'})

        report_params = {
            'message': '',
            'html_links': output_html_files,
            'direct_html_link_index': 0,
            'workspace_name': params['workspace_name'],
            'report_object_name': f'compare_annotations_{uuid.uuid4()}'}

        output = self.kbr.create_extended_report(report_params)

        return {'report_name': output['name'],
                'report_ref': output['ref']}

    def run(self, ctx, params):

        # collect some metadata
        self.get_genome(params['genome'])
        self.get_ontology_events(self.genome)
        self.get_translations()

        # make reports
        summary = self.summarize_gto(self.genome, self.translations)

        report = self.html_summary(params, summary)
        return report
