import os
import datetime
import logging
import json
import uuid
import yaml

from installed_clients.GenomeAnnotationAPIClient import GenomeAnnotationAPI
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.WorkspaceClient import Workspace as Workspace
from installed_clients.KBaseReportClient import KBaseReport

class CompareAnnotationsUtil:

    workdir     = 'tmp/work/'
    staging_dir = "/staging/"
    datadir     = "/kb/module/data/"

    def __init__(self, config):
        os.makedirs(self.workdir, exist_ok = True)
        self.config = config
        self.timestamp = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
        self.genes = {}
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.genome_api = GenomeAnnotationAPI(self.callback_url)
        self.dfu = DataFileUtil(self.callback_url)
        self.gfu = GenomeFileUtil(self.callback_url)
        self.kbr = KBaseReport(self.callback_url)
        self.ws_client = Workspace(config["workspace-url"])

    def run(self, ctx, params):
        logging.info("HELLO WORLD!!")
