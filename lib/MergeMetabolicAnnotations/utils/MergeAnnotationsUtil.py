import os
import datetime
import logging
import json
import uuid

from installed_clients.GenomeAnnotationAPIClient import GenomeAnnotationAPI
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.WorkspaceClient import Workspace as Workspace
from installed_clients.KBaseReportClient import KBaseReport

import MergeMetabolicAnnotations.utils.utils as mu


class MergeAnnotationsUtil:

    workdir = 'tmp/work/'
    staging_dir = "/staging/"
    datadir = "/kb/module/data/"

    translation_locations = {'ec': 'EBI_EC.ModelSEED.json',
                             'keggro': 'KEGG_RXN.ModelSEED.json',
                             'keggko': 'KEGG_KO.ModelSEED.json',
                             'metacyc': 'Metacyc_RXN.ModelSEED.json',
                             'SSO': 'SSO.ModelSEED.json',
                             'go': 'GO.ModelSEED.json'}
    ontology_lookup = {
        "ec": "EBI_EC_ontologyDictionary.json",
        "keggko": "KEGG_KO_ontologyDictionary.json",
        "keggro": "KEGG_RXN_ontologyDictionary.json",
        "metacyc": "MetaCyc_RXN_ontologyDictionary.json",
        "modelseed": "ModelSEED_RXN_ontologyDictionary.json",
        "go": "GO_ontologyDictionary.json",
        "sso": "SSO_ontologyDictionary.json"
    }

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
        self.weights = {}

    def get_ontology_events(self, genome_dict, params):
        if 'ontology_events' in genome_dict:
            for event, ontology in enumerate(genome_dict['ontology_events']):
                if ontology['description'] in params['annotations_to_merge'] or len(params['annotations_to_merge']) == 0:
                    self.events[event] = {}
                    for term in ontology:
                        self.events[event][term] = ontology[term]
                    if len(params['annotations_to_merge']) == 0:
                        self.weights[event] = 1  # If no annotations selected, weight them all equally
                    else:
                        for annot_item in params['annotations_to_merge']:
                            if annot_item['annotation_source'][0] == ontology['description']:
                                 self.weights[event] = annot_item['annotation_weight']
        else:
            logging.info("No ontology events in this genome!")

        logging.info(self.events)
        logging.info(self.weights)

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

    def merge_and_summarize_genome_dict(self, genome_dict, translations, params):
        # summary = {"genes": set(),
        #            "terms": set(),
        #            "rxns": set(),
        #            "orphan_terms": set()
        #            }
        summary = {}
        genes = {}

        summary[self.merged_ontology_event]= {"genes": set(),
                                              "terms": set(),
                                              "rxns": set(),
                                              "orphan_terms": set()
                                              }

        # add gene id to summary
        for feature in genome_dict['features']:
            gene_id = feature['id']
            reactions = {}

            # get ontology term
            if "ontology_terms" in feature:
                for type in feature['ontology_terms']:
                    term_dict = feature['ontology_terms'][type]

                    for term in term_dict:
                        for oe in term_dict[term]:

                            # is this ontology event in the user-selected list?
                            if oe in self.events:
                                if oe not in summary:
                                    summary[oe] = {"genes": set(),
                                                   "terms": set(),
                                                   "rxns": set(),
                                                   "orphan_terms": set()
                                                   }
                                summary[oe]["genes"].add(gene_id)

                                rxn = "none"

                                # get rxn
                                ontology_type = self.events[oe]['id']

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
                                    if term in genome_dict['ontologies_present']['SSO']:
                                        if genome_dict['ontologies_present']['SSO'][term] != 'Unknown':
                                            term = genome_dict['ontologies_present']['SSO'][term]

                                summary[oe]["terms"].add(term)

                                # convert terms to rxns
                                if term in translations[ontology_type]:
                                    rxn = translations[ontology_type][term]
                                    summary[oe]["rxns"].add(rxn)
                                else:
                                    summary[oe]["orphan_terms"].add(term)

                                # add up the weighted score for each reaction
                                if rxn in reactions:
                                    reactions[rxn] = reactions[rxn] + self.weights[oe]
                                else: reactions[rxn] = self.weights[oe]

                for rxn, score in reactions:
                    # keep rxn if score > threshold, and is max score OR user wants to keep all above threshold
                    if score >= params['annotation_threshold'] and (score == max(reactions.values()) or not params['keep_best_annotation_only']):
                        # collect gene_id:rxn in "genes" dict to be added using update_genome later
                        if gene_id not in genes:
                            genes[gene_id] = mu.Gene(gene_id)
                        genes[gene_id].addAnnotation(rxn)
                        genes[gene_id].valid=1
                        summary[self.merged_ontology_event]["genes"].add(gene_id)
                        summary[self.merged_ontology_event]["terms"].add(rxn)
                        summary[self.merged_ontology_event]["rxns"].add(rxn)

        # # caution: json.dump cannot handle objects wit sets in them?
        # with open(os.path.join(self.scratch, "summary_dump.json"), 'w') as outfile:
        #     json.dump(summary, outfile, indent=2)

        return genes, summary

    # def merge_and_summarize_genome_dict(self, genome_dict, translations, params):
    #     summary = {"genes": {},
    #                "terms": {},
    #                "rxns": {},
    #                "ontology_events": {},
    #                "orphan_terms": {}
    #                }
    #     genes = {}
    #
    #     # add ontology events
    #     for oe in self.events:
    #         summary['ontology_events'][oe] = self.events[oe]
    #
    #     # add gene id to summary
    #     for feature in genome_dict['features']:
    #         gene_id = feature['id']
    #         summary["genes"][gene_id] = {"terms": {},
    #                                      "rxns": {}
    #                                      }
    #         reactions = {}
    #
    #         # get ontology term
    #         if "ontology_terms" in feature:
    #             for type in feature['ontology_terms']:
    #                 term_dict = feature['ontology_terms'][type]
    #
    #                 for term in term_dict:
    #                     for oe in term_dict[term]:
    #
    #                         # is this ontology event in the user-selected list?
    #                         if oe in self.events:
    #
    #                             rxn = "none"
    #
    #                             # get rxn
    #                             ontology_type = summary['ontology_events'][oe]['id']
    #
    #                             # fix metacyc terms
    #                             if ontology_type == 'metacyc':
    #                                 if term.startswith("META:"):
    #                                     term = term.replace('META:', '')
    #
    #                             # fix go terms
    #                             if ontology_type == 'go':
    #                                 if term.startswith("GO:"):
    #                                     term = term.replace('GO:', '')
    #
    #                             # fix SSO terms
    #                             if ontology_type == 'SSO':
    #                                 if term in genome_dict['ontologies_present']['SSO']:
    #                                     if genome_dict['ontologies_present']['SSO'][term] != 'Unknown':
    #                                         term = genome_dict['ontologies_present']['SSO'][term]
    #
    #                             # convert terms to rxns
    #                             if term in translations[ontology_type]:
    #                                 rxn = translations[ontology_type][term]
    #                             else:
    #                                 if oe in summary["orphan_terms"]:
    #                                     summary["orphan_terms"][oe].append(term)
    #                                     summary["orphan_terms"][oe] = list(
    #                                         set(summary["orphan_terms"][oe]))
    #                                 else:
    #                                     summary["orphan_terms"][oe] = [term]
    #
    #                             # terms
    #                             if term in summary["genes"][gene_id]['terms']:
    #                                 summary["genes"][gene_id]['terms'][term].append(oe)
    #                             else:
    #                                 summary["genes"][gene_id]['terms'][term] = [oe]
    #
    #                             if term in summary['terms']:
    #                                 summary['terms'][term].append(oe)
    #                                 summary['terms'][term] = list(set(summary['terms'][term]))
    #                             else:
    #                                 summary['terms'][term] = [oe]
    #
    #                             # rxns
    #                             if rxn != "none":
    #                                 if rxn in summary["genes"][gene_id]['rxns']:
    #                                     summary["genes"][gene_id]['rxns'][rxn].append(oe)
    #                                 else:
    #                                     summary["genes"][gene_id]['rxns'][rxn] = [oe]
    #
    #                                 if rxn in summary['rxns']:
    #                                     summary['rxns'][rxn].append(oe)
    #                                     summary['rxns'][rxn] = list(set(summary['rxns'][rxn]))
    #                                 else:
    #                                     summary['rxns'][rxn] = [oe]
    #
    #                                 # add up the weighted score for each reaction
    #                                 if rxn in reactions:
    #                                     reactions[rxn] = reactions[rxn] + self.weights[oe]
    #                                 else: reactions[rxn] = self.weights[oe]
    #
    #             for rxn, score in reactions:
    #                 # keep rxn if score > threshold, and is max score OR user wants to keep all above threshold
    #                 if score >= params['annotation_threshold'] and (score == max(reactions.values()) or not params['keep_best_annotation_only']):
    #                     # collect gene_id:rxn in "genes" dict to be added using update_genome later
    #                     if gene_id not in genes:
    #                         genes[gene_id] = mu.Gene(gene_id)
    #                     genes[gene_id].addAnnotation(rxn)
    #
    #     with open(os.path.join(self.scratch, "summary_dump.json"), 'w') as outfile:
    #         json.dump(summary, outfile, indent=2)
    #
    #     return genes, summary

    def html_summary(self, params, summary):

        # # convert genome_dict summary for this report
        # html_summary_report = {}
        #
        # for oe in summary:
        #     html_summary_report[oe] = {"gene": [], "term": [], "rxn": []}
        #
        # for gene in summary["genes"]:
        #     for term in summary["genes"][gene]['terms']:
        #         for oe in summary["genes"][gene]['terms'][term]:
        #             html_summary_report[oe]['gene'].append(gene)
        #             html_summary_report[oe]['term'].append(term)
        #
        #             html_summary_report[oe]['gene'] = list(set(html_summary_report[oe]['gene']))
        #             html_summary_report[oe]['term'] = list(set(html_summary_report[oe]['term']))
        #
        #     for rxn in summary["genes"][gene]['rxns']:
        #         for oe in summary["genes"][gene]['rxns'][rxn]:
        #             html_summary_report[oe]['rxn'].append(rxn)
        #             html_summary_report[oe]['gene'].append(gene)
        #
        #             html_summary_report[oe]['rxn'] = list(set(html_summary_report[oe]['rxn']))
        #             html_summary_report[oe]['gene'] = list(set(html_summary_report[oe]['gene']))

        output_html_files = list()

        # Make report directory and copy over files
        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        os.mkdir(output_directory)
        result_file_path = os.path.join(output_directory, 'merge_annotations_summary.html')

        # make html
        table_lines = []
        table_lines.append(f'<h2>Merge Annotations</h2>')
        table_lines.append(f'<h3>Summary</h3>')
        table_lines.append(
            '<table cellspacing="0" cellpadding="3" border="1"><tr><th>EVENT</th><th>DESCRIPTION</th><th>TYPE</th><th>GENES</th><th>TERMS</th><th>RXNS</th></tr>')
        for event in sorted(summary.keys()):
            # RAST/PROKKA don't have descriptions, but they have methods
            description = self.events[event].get('description', self.events[event]['method'])
            type = self.events[event]['id']
            table_lines.append('<tr><td>' + str(event) + \
                              '</td><td>' + description + \
                              '</td><td>' + type + \
                              '</td><td>' + str(len(summary[event]["genes"])) + \
                              '</td><td>' + str(len(summary[event]["terms"])) + \
                              '</td><td>' + str(len(summary[event]["rxns"])) + \
                              '</td></tr>')
        table_lines.append('</table>')

        # Write to file
        with open(result_file_path, 'w') as result_file:
            for line in table_lines:
                result_file.write(line + "\n")

        output_html_files.append(
            {'path': output_directory,
             'name': os.path.basename(result_file_path),
             'description': 'Summary Report'})

        # # bokeh plots
        # totals_file_path = os.path.join(output_directory, 'totals.html')
        # output_file(totals_file_path, title="Totals")
        # totals = self.plot_totals(summary)
        # save(totals)
        # output_html_files.append(
        #     {'path': output_directory,
        #      'name': os.path.basename(totals_file_path),
        #      'description': 'Ontology Totals'})
        #
        # csc_file_path = os.path.join(output_directory, 'csc.html')
        # output_file(csc_file_path, title="CSC")
        # csc = self.plot_csc2(summary)
        # save(csc)
        # output_html_files.append(
        #     {'path': output_directory,
        #      'name': os.path.basename(csc_file_path),
        #      'description': 'Cumulative Sum Plot'})

        # finalize html reports
        report_params = {
            'message': '',
            'html_links': output_html_files,
            'direct_html_link_index': 0,
            'workspace_name': params['workspace_name'],
            'report_object_name': f'merge_annotations_{uuid.uuid4()}'}

        output = self.kbr.create_extended_report(report_params)

        return {'report_name': output['name'],
                'report_ref': output['ref']}

    def run(self, ctx, params):

        # read and prepare objects/files
        mu.validate()
        self.sso_ref = mu.get_sso_data("KBaseOntology/seed_subsystem_ontology", \
                       self.ws_client)

        # get genome
        get_genome_results = mu.get_genome(params['genome'], self.genome_api)
        genome_dict = get_genome_results[0]
        self.genome_full = get_genome_results[1]

        ontology_dict = mu.get_ontology_dict(params['ontology'], self.datadir, self.ontology_lookup)

        # collect some metadata
        self.get_ontology_events(genome_dict, params)
        self.get_translations()

        genome_dict = mu.add_ontology_event(genome_dict,
                                            {"ontology": "modelseed",
                                             "description": "Merged annotations"},
                                            self.sso_ref, self.timestamp,
                                            "Merge Annotations")
        self.merged_ontology_event = len(genome_dict['ontology_events']) - 1

        # make reports
        genes, summary = self.merge_and_summarize_genome_dict(genome_dict, \
                         self.translations, params)

#        report = self.html_summary(params, summary)

        # # process - do we need this? If so, also need to set ontology_dict
        # for gene in self.genes:
        #     self.genes[gene].validateGeneID(genome_dict)
        #     # logging.info(self.genes[gene].id)
        for gene in self.genes:
            self.genes[gene].validateAnnotationID(ontology_dict, "modelseed"")

        genome_dict = mu.update_genome(genome_dict, "modelseed", \
                      genes, self.merged_ontology_event)

        # overwrite object with new genome
        self.genome_full['data'] = genome_dict

        prov = ctx.provenance()
        info = self.gfu.save_one_genome({'workspace': params['workspace_name'],
                                         'name': params['output_name'],
                                         'data': self.genome_full['data'],
                                         'provenance': prov})['info']

        genome_ref = str(info[6]) + '/' + str(info[0]) + '/' + str(info[4])
        logging.info("*** Genome ID: " + str(genome_ref))

        report = self.generate_report(params, genome_ref)

        return report
