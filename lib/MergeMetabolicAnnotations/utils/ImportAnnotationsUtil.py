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

    def run(self, ctx, params):

        ontology = f.df_to_ontology(params)

        add_ontology_results = self.anno_api.add_annotation_ontology_events({
            "input_ref": params['genome'],
            "output_name": params['output_name'],
            "input_workspace": params['workspace_name'],
            "workspace-url": self.config["workspace-url"],
            "events": [ontology],
            "timestamp": self.timestamp,
            "output_workspace": params['workspace_name'],
            "save": 1
        })

        # get the new list of events to make a table
        get_ontology_results = self.anno_api.get_annotation_ontology_events({
            "input_ref": add_ontology_results['output_ref'],
            "workspace-url": self.config["workspace-url"]
        })

        ontology_selected = f.filter_selected_ontologies(
            get_ontology_results, params, workflow="unique")
        with open(os.path.join(self.scratch, "get_ontology_dump.json"), 'w') as outfile:
            json.dump(ontology_selected, outfile, indent=2)

        # make report
        html_reports = []
        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        os.mkdir(output_directory)

        html_reports.append(f.html_add_ontology_summary(
            params, ontology, add_ontology_results, output_directory))
        event_summary = f.get_event_lists(ontology_selected)
        html_reports = f.compare_report_stack(html_reports, event_summary, output_directory)

        # finalize html reports
        report_params = {
            'message': '',
            'html_links': html_reports,
            'direct_html_link_index': 0,
            'objects_created': [{'ref': add_ontology_results["output_ref"], 'description': 'Genome with imported annotations'}],
            'workspace_name': params['workspace_name'],
            'report_object_name': f'import_annotations_{uuid.uuid4()}'}

        report_output = self.kbr.create_extended_report(report_params)

        return {'output_genome_ref': add_ontology_results["output_ref"],
                'report_name': report_output['name'],
                'report_ref': report_output['ref']}
