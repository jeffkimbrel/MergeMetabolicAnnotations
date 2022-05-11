import os
import datetime
import logging
import json
import uuid

from installed_clients.WorkspaceClient import Workspace as Workspace
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.annotation_ontology_apiServiceClient import annotation_ontology_api

import MergeMetabolicAnnotations.utils.functions as f


class CompareAnnotationsUtil:

    def __init__(self, config):
        self.config = config
        self.timestamp = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.kbr = KBaseReport(self.callback_url)
        self.anno_api = annotation_ontology_api()
        self.ws_client = Workspace(config["workspace-url"])

    def run(self, ctx, params):
        get_ontology_results = self.anno_api.get_annotation_ontology_events({
            "input_ref": params['genome'],
            "workspace-url": self.config["workspace-url"]
        })

        ontology_selected = f.filter_selected_ontologies(
            get_ontology_results, params, workflow="compare")
        
        # make reports
        html_reports = []
        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        os.mkdir(output_directory)

        event_summary = f.get_event_lists(ontology_selected)
        html_reports = f.compare_report_stack(html_reports, event_summary, output_directory)

        # add json dump
        with open(os.path.join(output_directory, "get_ontology_dump.json"), 'w') as outfile:
            json.dump(ontology_selected, outfile, indent=2)
            
        html_reports.append({'path': output_directory,
            'name': "get_ontology_dump.json",
            'description': 'all results'})

        # finalize html reports
        report_params = {
            'message': '',
            'html_links': html_reports,
            'direct_html_link_index': 0,
            'workspace_name': params['workspace_name'],
            'report_object_name': f'compare_annotations_{uuid.uuid4()}'}

        report_output = self.kbr.create_extended_report(report_params)

        return {'report_name': report_output['name'],
                'report_ref': report_output['ref']}
