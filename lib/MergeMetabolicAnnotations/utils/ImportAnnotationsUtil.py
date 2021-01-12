# import sys
import os
import datetime
import logging
import json
import uuid

from installed_clients.WorkspaceClient import Workspace as Workspace
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.annotation_ontology_apiServiceClient import annotation_ontology_api

import MergeMetabolicAnnotations.utils.functions as f


class ImportAnnotationsUtil:

    def __init__(self, config):
        self.config = config
        self.timestamp = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.ws_client = Workspace(config["workspace-url"])
        self.anno_api = annotation_ontology_api()
        self.kbr = KBaseReport(self.callback_url)

    def generate_report(self, params, ontology, output):

        # Make report directory and copy over files
        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        os.mkdir(output_directory)
        result_file_path = os.path.join(output_directory, 'import_annotations_summary.html')

        report = []
        report.append(f'<h3>Import Annotations Summary</h3>')

        report.append(f'Import App version: {ontology["method_version"]}<br>')
        report.append(f'Timestamp: {ontology["timestamp"]}<br>')

        report.append(f'Annotations file: {params["annotation_file"]}<br>')
        report.append(f'Ontology ID: {params["ontology"]}<br>')
        report.append(f'Description: {params["description"]}<br>')

        report.append(f'Input Ref: {params["genome"]}<br>')
        report.append(f'Output Ref: {output["output_ref"]}<br><br>')

        report.append(f'Features in annotations file: {len(ontology["ontology_terms"])}<br>')
        report.append(f'Features (found): {output["ftrs_found"]}<br>')
        report.append(f'Features (not found): {len(output["ftrs_not_found"])}<br>')

        if len(output["ftrs_not_found"]) > 0:
            report.append(
                f'These genes were not found in the genome: <br>{(", ").join(output["ftrs_not_found"])}<br>')

        # Write to file
        with open(result_file_path, 'w') as result_file:
            for line in report:
                result_file.write(line + "\n")

        output_html_files = [
            {'path': output_directory,
             'name': os.path.basename(result_file_path),
             'description': 'HTML report for import_annotations app'}
        ]

        report_params = {
            'message': '',
            'html_links': output_html_files,
            'direct_html_link_index': 0,
            'objects_created': [{'ref': output["output_ref"], 'description': 'Genome with imported annotations'}],
            'workspace_name': params['workspace_name'],
            'report_object_name': f'import_annotations_{uuid.uuid4()}'}

        report_output = self.kbr.create_extended_report(report_params)

        return {'output_genome_ref': output["output_ref"],
                'report_name': report_output['name'],
                'report_ref': report_output['ref']}

    def run(self, ctx, params):

        ontology = f.df_to_ontology(params)

        output = self.anno_api.add_annotation_ontology_events({
            "input_ref": params['genome'],
            "output_name": params['output_name'],
            "input_workspace": params['workspace_name'],
            "workspace-url": self.config["workspace-url"],
            "events": [ontology],
            "timestamp": self.timestamp,
            "output_workspace": params['workspace_name'],
            "save": 1
        })

        report = self.generate_report(params, ontology, output)

        return report
