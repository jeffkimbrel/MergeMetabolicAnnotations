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

import MergeMetabolicAnnotations.utils.utils as mu


class ImportBulkAnnotationsUtil:

    workdir = 'tmp/work/'
    staging_dir = "/staging/"
    datadir = "/kb/module/data/"

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

        get_genome_results = mu.get_genome(params['genome'], self.genome_api)
        genome_dict = get_genome_results[0]
        self.genome_full = get_genome_results[1]

        bulk_annotations = mu.get_bulk_annotations_file(params, self.staging_dir)

        mu.validate_bulk(bulk_annotations)

        # identify and process each pair of descriptions/ontologies
        pairs = mu.get_description_ontology_pairs(bulk_annotations)
        for index, row in pairs.iterrows():
            # logging.info("---")
            #logging.info(row['description'] + "\t" + row['ontology'])

            pair_params = params
            pair_params['description'] = row['description']
            pair_params['ontology'] = row['ontology']

            genes = {}
            sso_ref = mu.get_sso_data(pair_params['ontology'], self.ws_client)
            genome_dict = mu.add_ontology_event(
                genome_dict, pair_params, sso_ref, self.timestamp, "Import Annotations")
            current_ontology_event = len(genome_dict['ontology_events']) - 1
            ontology_dict = mu.get_ontology_dict(
                pair_params['ontology'], self.datadir, self.ontology_lookup)

            annotations = bulk_annotations[(bulk_annotations.description == row['description'])
                                           & (bulk_annotations.ontology == row['ontology'])][['gene', 'term']]

            genes = mu.annotations_to_genes(annotations, genes)

            for gene in genes:
                genes[gene].validateGeneID(genome_dict)
                genes[gene].validateAnnotationID(ontology_dict, pair_params['ontology'])

            genome_dict = mu.update_genome(
                genome_dict, pair_params['ontology'], genes, current_ontology_event)

        # overwrite object with new genome
        self.genome_full['data'] = genome_dict

        prov = ctx.provenance()
        info = self.gfu.save_one_genome({'workspace': params['workspace_name'],
                                         'name': params['output_name'],
                                         'data': self.genome_full['data'],
                                         'provenance': prov})['info']

        genome_ref = str(info[6]) + '/' + str(info[0]) + '/' + str(info[4])
        logging.info("*** Genome ID: " + str(genome_ref))

        # report = self.generate_report(params, genome_ref)
        #
        # return report

    def run_old(self, ctx, params):
        """
        The main run called by the implementation file.
        """

        # read and prepare objects/files
        mu.validate()
        self.sso_ref = mu.get_sso_data(params['ontology'], self.ws_client)

        # get genome
        get_genome_results = mu.get_genome(params['genome'], self.genome_api)
        genome_dict = get_genome_results[0]
        self.genome_full = get_genome_results[1]

        ontology_dict = mu.get_ontology_dict(params['ontology'], self.datadir, self.ontology_lookup)
        annotations = mu.get_annotations_file(params, self.staging_dir)

        self.genes = mu.annotations_to_genes(annotations, self.genes)
        genome_dict = mu.add_ontology_event(genome_dict, params, self.sso_ref, self.timestamp)
        self.current_ontology_event = len(genome_dict['ontology_events']) - 1

        # process
        for gene in self.genes:
            self.genes[gene].validateGeneID(genome_dict)
            self.genes[gene].validateAnnotationID(ontology_dict, params['ontology'])

        genome_dict = mu.update_genome(
            genome_dict, params['ontology'], self.genes, self.current_ontology_event)

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
