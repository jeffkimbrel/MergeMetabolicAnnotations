import os
import datetime
import logging
import json
import uuid
import pandas as pd
from collections import Counter

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
        self.genes = {}

    def get_ontology_events(self, params):
        if 'ontology_events' in self.genome:

            for event, ontology in enumerate(self.genome['ontology_events']):

                # fix some legacy problems
                if 'description' not in ontology:
                    ontology['description'] = ontology['method']
                ontology["id"] = mu.legacy_fix(ontology["id"])

                if len(params['annotations_to_merge']) == 0:
                    self.weights[event] = 1
                    self.events[event] = {}
                    for term in ontology:
                        self.events[event][term] = ontology[term]

                else:
                    for annotations_to_merge in params['annotations_to_merge']:
                        if ontology['description'] in annotations_to_merge['annotation_source']:
                            self.events[event] = {}
                            self.weights[event] = annotations_to_merge['annotation_weight']
                            for term in ontology:
                                self.events[event][term] = ontology[term]

        else:
            logging.info("No ontology events in this genome!")

    def merge_annotations(self):
        merged_annotations = {}

        # add gene id to summary
        for feature in self.genome['features']:
            gene_id = feature['id']
            merged_annotations[gene_id] = {}

            # get ontology term
            if "ontology_terms" in feature:
                for type in feature['ontology_terms']:
                    term_dict = feature['ontology_terms'][type]

                    # fix potential legacy problems after getting feature
                    type = mu.legacy_fix(type)

                    for term in term_dict:
                        # logging.info(term)
                        # logging.info(mu.standardize_annotation(term, type))

                        for ontology_event in term_dict[term]:

                            # is this ontology event in the user-selected list?
                            if ontology_event in self.events:

                                rxn = "none"

                                # convert terms to rxns
                                standardized_term = mu.standardize_annotation(term, type)

                                if standardized_term in self.translations[type]:
                                    rxn = self.translations[type][standardized_term]

                                if rxn != "none":
                                    if rxn in merged_annotations[gene_id]:
                                        merged_annotations[gene_id][rxn]['events'].append(
                                            ontology_event)

                                        # clean up duplicates... eg old versions of prokka added many of the same reaction
                                        merged_annotations[gene_id][rxn]['events'] = list(
                                            set(merged_annotations[gene_id][rxn]['events']))
                                    else:
                                        merged_annotations[gene_id][rxn] = {'events': []}
                                        merged_annotations[gene_id][rxn]['events'] = [
                                            ontology_event]

        return merged_annotations

    def score_annotations(self, annotations, threshold, best_only):
        '''
        returns a pandas dataframe suitable for import annotations
        '''

        df = pd.DataFrame(columns=['gene', 'term', 'events', 'score'])

        for gene_id in annotations:

            # get total score of each rxn, save to 'score_total'
            for rxn in annotations[gene_id]:
                annotations[gene_id][rxn]['score_total'] = 0
                for ontology_event in annotations[gene_id][rxn]['events']:
                    annotations[gene_id][rxn]['score_total'] += self.weights[ontology_event]

            # get list of the best rxn or rxns
            '''
            best_score = max()

            for rxn in annotations[gene_id]:

                if annotations[gene_id][rxn]['score_total'] >= threshold:
                    if best_only == "all" or rxn in best_rxn:
                        annotations[gene_id][rxn]['passed'] = 1
                        df = df.append(
                            pd.Series(data={'gene': gene_id, 'term': rxn, 'events': annotations[gene_id][rxn]['events'], 'score': annotations[gene_id][rxn]['score_total']}), ignore_index=True)

                    else:
                        annotations[gene_id][rxn]['passed'] = 0
                        # df = df.append(
                        #     pd.Series(data={'gene': gene_id, 'term': rxn, 'events': annotations[gene_id][rxn]['events'], 'score': annotations[gene_id][rxn]['score_total']}), ignore_index=True)

            #        annotations[gene_id][rxn]['passed'] = 0
            '''

        # with open(os.path.join(self.scratch, "scored.json"), 'w') as outfile:
        #     json.dump(annotations, outfile, indent=2)

        df.to_csv(os.path.join(self.scratch, "scored.txt"), sep="\t", index=False)

        return df

    def html_summary(self, params, summary):

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
            table_lines.append('<tr><td>' + str(event) +
                               '</td><td>' + description +
                               '</td><td>' + type +
                               '</td><td>' + str(len(summary[event]["genes"])) +
                               '</td><td>' + str(len(summary[event]["terms"])) +
                               '</td><td>' + str(len(summary[event]["rxns"])) +
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

    def generate_report(self, params, genome_ref):
        """
        Reads in the results from the summary method, and creates the html
        report.

        This is just a copy/paste of the report from the import app
        """

        summary = mu.summarize(params, self.genes)

        output_html_files = list()

        # Make report directory and copy over files
        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        os.mkdir(output_directory)
        result_file_path = os.path.join(output_directory, 'import_annotations_summary.html')

        # Build HTML tables for results
        table_lines = []
        table_lines.append(f'<h2>Import Annotations</h2>')
        table_lines.append(f'<h3>Summary</h3>')
        table_lines.append(
            '<table cellspacing="0" cellpadding="3" border="1"><tr><th>TYPE</th><th>VALID</th><th>INVALID</th></tr>')
        table_lines.append('<tr><td>GENES</td><td>' + str(
            len(summary['valid_genes'])) + '</td><td>' + str(len(summary['invalid_genes'])) + '</td></tr>')
        table_lines.append('<tr><td>TERMS</td><td>' + str(
            len(summary['valid_terms'])) + '</td><td>' + str(len(summary['invalid_terms'])) + '</td></tr>')
        table_lines.append('</table>')

        if len(summary['invalid_genes']) > 0:
            table_lines.append(f'<h3>Invalid Genes</h3>')
            table_lines.append(
                '<i>These are locus_tags not identified in the genome object. Frequency shown in parentheses.</i><br><br>')

            invalid_genes_count = dict(Counter(summary['invalid_genes']))

            for gene in sorted(invalid_genes_count.keys()):
                gene_count = gene + '\t(' + str(invalid_genes_count[gene]) + ')'
                table_lines.append(gene_count + '<br>')

        if len(summary['invalid_terms']) > 0:
            table_lines.append(f'<h3>Invalid Terms</h3>')
            table_lines.append(
                '<i>These are ontology terms not found in the ontology dictionary. Frequency shown in parentheses.</i><br><br>')

            invalid_terms_count = dict(Counter(summary['invalid_terms']))

            for term in sorted(invalid_terms_count.keys()):
                term_count = term + '\t(' + str(invalid_terms_count[term]) + ')'
                table_lines.append(term_count + '<br>')

        # Write to file
        with open(result_file_path, 'w') as result_file:
            for line in table_lines:
                result_file.write(line + "\n")

        output_html_files.append(
            {'path': output_directory,
             'name': os.path.basename(result_file_path),
             'description': 'HTML report for import_annotations app'})

        report_params = {
            'message': '',
            'html_links': output_html_files,
            'direct_html_link_index': 0,
            'objects_created': [{'ref': genome_ref, 'description': 'Genome with imported annotations'}],
            'workspace_name': params['workspace_name'],
            'report_object_name': f'import_annotations_{uuid.uuid4()}'}

        output = self.kbr.create_extended_report(report_params)

        return {'output_genome_ref': genome_ref,
                'report_name': output['name'],
                'report_ref': output['ref']}

    def run(self, ctx, params):
        params['ontology'] = 'MSRXN'  # just in case it doesn't get set

        self.genome = mu.get_genome(params['genome'], self.genome_api)

        self.get_ontology_events(params)
        self.translations = mu.get_translations(self.datadir)

        merged_annotations = self.merge_annotations()
        scored_annotations = self.score_annotations(
            merged_annotations,
            params['annotation_threshold'],
            params['keep_best_annotation_only'])

        ontology_dict = mu.get_ontology_dict('MSRXN',
                                             self.datadir,
                                             mu.ontology_lookup)

        # get list of uploaded annotation terms
        annotations = mu.get_annotations_file(params, self.staging_dir, pass_df=scored_annotations)
        self.genes = mu.annotations_to_genes(annotations, self.genes)

        self.genome = mu.add_ontology_event(
            self.genome, params, self.timestamp, "Merge Annotations")

        # fix missing descriptions
        o_counter = 0
        for ontology_event in self.genome['ontology_events']:
            if 'description' not in ontology_event:
                ontology_event['description'] = ontology_event['method']
            self.genome['ontology_events'][o_counter] = ontology_event
            o_counter += 1

        self.current_ontology_event = len(self.genome['ontology_events']) - 1

        # process
        for gene in self.genes:
            self.genes[gene].validate_gene_ID(self.genome)
            self.genes[gene].validate_annotation_ID(ontology_dict, 'MSRXN')

        self.genome = mu.update_genome(
            self.genome,
            'MSRXN',
            self.genes,
            self.current_ontology_event)

        info = self.gfu.save_one_genome({'workspace': params['workspace_name'],
                                         'name': params['output_name'],
                                         'data': self.genome,
                                         'provenance': ctx.provenance()})['info']

        genome_ref = str(info[6]) + '/' + str(info[0]) + '/' + str(info[4])
        logging.info("*** Genome ID: " + str(genome_ref))

        report = self.generate_report(params, genome_ref)

        return report

    # def run_old(self, ctx, params):
    #
    #     # read and prepare objects/files
    #     mu.validate()
    #
    #     # get genome
    #     self.genome = mu.get_genome(params['genome'], self.genome_api)
    #
    #     ontology_dict = mu.get_ontology_dict(params['ontology'], self.datadir, mu.ontology_lookup)
    #
    #     # collect some metadata
    #     self.get_ontology_events(params)
    #     logging.info(str(self.events))
    #     # self.get_translations()
    #
    #     self.genome = mu.add_ontology_event(self.genome,
    #                                         {"ontology": "MSRXN",
    #                                          "description": "Merged annotations"},
    #                                         self.timestamp,
    #                                         "Merge Annotations")
    #
    #     self.merged_ontology_event = len(self.genome['ontology_events']) - 1
    #
    #     # make reports
    #     genes, summary = self.merge_and_summarize_genome_dict(params)
    #
    #     for gene in genes:
    #         genes[gene].validateAnnotationID(ontology_dict, "MSRXN")
    #
    #     self.genome = mu.update_genome(self.genome, "MSRXN",
    #                                    genes, self.merged_ontology_event)
    #
    #     info = self.gfu.save_one_genome({'workspace': params['workspace_name'],
    #                                      'name': params['output_name'],
    #                                      'data': self.genome,
    #                                      'provenance': ctx.provenance()})['info']
    #
    #     genome_ref = str(info[6]) + '/' + str(info[0]) + '/' + str(info[4])
    #     logging.info("*** Genome ID: " + str(genome_ref))
    #
    #     # add event for reporting
    #     self.events[self.merged_ontology_event] = {"ontology": "MSRXN",
    #                                                "description": "Merged annotations",
    #                                                "method": "Merged annotations",
    #                                                "id": "MSRXN"}
    #     report = self.html_summary(params, summary)
    #
    #     return report
