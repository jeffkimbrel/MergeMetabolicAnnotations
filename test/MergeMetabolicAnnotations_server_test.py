# -*- coding: utf-8 -*-
import os
import time
import unittest
from configparser import ConfigParser

from MergeMetabolicAnnotations.MergeMetabolicAnnotationsImpl import MergeMetabolicAnnotations
from MergeMetabolicAnnotations.MergeMetabolicAnnotationsServer import MethodContext
from MergeMetabolicAnnotations.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace


class MergeMetabolicAnnotationsTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('MergeMetabolicAnnotations'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'MergeMetabolicAnnotations',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = MergeMetabolicAnnotations(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_ContigFilter_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    def test_your_method(self):
        # Prepare test objects in workspace if needed using
        # self.getWsClient().save_objects({'workspace': self.getWsName(),
        #                                  'objects': []})
        #
        # Run your method by
        # ret = self.getImpl().your_method(self.getContext(), parameters...)
        #
        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods

        # # import app
        # params_import = {
        #     "debug": True,
        #     "ontology": "go",
        #     "annotation_file": "PT19DW_GO.txt",
        #     "description": "test",
        #     "genome": "30554/14/1",
        #     "output_name": "testGenome_2",
        #     "workspace_name": self.wsName
        # }
        # ret = self.serviceImpl.import_annotations(self.ctx, params_import)

        # compare app
        params_compare = {
            "debug": True,
            "genome": "30554/17/1",
            "output_name": "compareGenome_temp",
            "annotations_to_compare": [],
            "workspace_name": self.wsName
        }
        ret = self.serviceImpl.compare_metabolic_annotations(self.ctx, params_compare)
