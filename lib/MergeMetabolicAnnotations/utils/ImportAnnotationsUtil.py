import sys
import os
import datetime
import logging
import json
import uuid
from collections import Counter

from installed_clients.GenomeAnnotationAPIClient import GenomeAnnotationAPI
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.WorkspaceClient import Workspace as Workspace
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.annotation_ontology_apiServiceClient import annotation_ontology_api

import MergeMetabolicAnnotations.utils.utils as mu
import MergeMetabolicAnnotations.utils.functions as f

#####


class ImportAnnotationsUtil:

    workdir = 'tmp/work/'
    staging_dir = "/staging/"
    datadir = "/kb/module/data/"

    def __init__(self, config):
        os.makedirs(self.workdir, exist_ok=True)
        self.config = config
        self.timestamp = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.workspace_url = config["workspace-url"]
        self.ws_client = Workspace(config["workspace-url"])
        self.anno_api = annotation_ontology_api()

        # self.genome_api = GenomeAnnotationAPI(self.callback_url)
        # self.dfu = DataFileUtil(self.callback_url)
        # self.gfu = GenomeFileUtil(self.callback_url)
        self.kbr = KBaseReport(self.callback_url)

        self.genes = {}

    def generate_report(self, params, genome_ref):
        """
        Reads in the results from the summary method, and creates the html
        report.
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
        '''
        The main run called by the implementation file.
        '''

        ontology = f.df_to_ontology(params, self.staging_dir)

        with open(os.path.join(self.scratch, "new_ontology_API_dump.json"), 'w') as outfile:
            json.dump(ontology, outfile, indent=2)

        output = self.anno_api.add_annotation_ontology_events({
            "input_ref": params['genome'],
            "output_name": params['output_name'],
            "input_workspace": self.workspace_url + params['workspace_name'],
            "events": [ontology],
            "timestamp": self.timestamp,
            "output_workspace": self.workspace_url + params['workspace_name'],
            "save": 1
        })

        logging.info(str(output))

        # DYNAMIC SERVICE
        #
        # logging.info(self.anno_api)
        # ontology = self.anno_api.get_annotation_ontology_events({
        #     "input_ref": params['genome']
        # })
        # with open(os.path.join(self.scratch, "ontology_API_dump.json"), 'w') as outfile:
        #     json.dump(ontology, outfile, indent=2)
        #
        # sys.exit()
        # #####
        #
        # # read and prepare imported file
        # mu.validate()
        #
        # # get genome object, store as a dictionary
        # self.genome = mu.get_genome(params['genome'], self.genome_api)
        #
        # # get ontology dictionary
        # ontology_dict = mu.get_ontology_dict(params['ontology'],
        #                                      self.datadir,
        #                                      mu.ontology_lookup)
        #
        # # get list of uploaded annotation terms
        # annotations = mu.get_annotations_file(params, self.staging_dir)
        #
        # # add annotation terms to gene class
        # self.genes = mu.annotations_to_genes(annotations, self.genes)
        #
        # self.genome = mu.add_ontology_event(
        #     self.genome, params, self.timestamp, "Import Annotations")
        #
        # # fix missing descriptions
        # o_counter = 0
        # for ontology_event in self.genome['ontology_events']:
        #     if 'description' not in ontology_event:
        #         ontology_event['description'] = ontology_event['method']
        #     self.genome['ontology_events'][o_counter] = ontology_event
        #     o_counter += 1
        #
        # self.current_ontology_event = len(self.genome['ontology_events']) - 1
        #
        # # process
        # for gene in self.genes:
        #     self.genes[gene].validate_gene_ID(self.genome)
        #     self.genes[gene].validate_annotation_ID(ontology_dict, params['ontology'])
        #
        # self.genome = mu.update_genome(
        #     self.genome,
        #     params['ontology'],
        #     self.genes,
        #     self.current_ontology_event)
        #
        # info = self.gfu.save_one_genome({'workspace': params['workspace_name'],
        #                                  'name': params['output_name'],
        #                                  'data': self.genome,
        #                                  'provenance': ctx.provenance()})['info']
        #
        # genome_ref = str(info[6]) + '/' + str(info[0]) + '/' + str(info[4])
        # logging.info("*** Genome ID: " + str(genome_ref))
        #
        # report = self.generate_report(params, genome_ref)
        #
        # return report
