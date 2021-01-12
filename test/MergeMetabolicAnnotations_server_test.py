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
        cls.wsName = "test_merge_" + str(suffix)
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

        # import app
        params_import = {
            "debug": True,
            "ontology": "KO",
            "annotation_file": "K_algicida_OT-1_protein_IDs.faa.kofam93.txt",
            "description": "annotation_API_test",
            "genome": "44643/3/1",
            "output_name": "import_genome",
            "workspace_name": self.wsName
        }
        ret = self.serviceImpl.import_annotations(self.ctx, params_import)

        # bulk import app
        # params_import = {
        #     "debug": True,
        #     "annotation_file": "pt32_bulk.txt",
        #     "description": "test",
        #     "genome": "30128/3/1",
        #     "output_name": "import_bulk_genome",
        #     "workspace_name": self.wsName
        # }
        # ret = self.serviceImpl.import_bulk_annotations(self.ctx, params_import)

        # compare app
        # params_compare = {
        #     "debug": True,
        #     "genome": "23001/49/1",
        #     "output_name": "compareGenome_temp",
        #     "annotations_to_compare": [],
        #     # "annotations_to_compare": ["ProkkaAnnotation:2.1.5:EC:2020_10_13_20_07_38", "KOFAM93", "annotate_genome"],
        #     "workspace_name": self.wsName
        # }
        # ret = self.serviceImpl.compare_metabolic_annotations(self.ctx, params_compare)

        # merge app
        # params_merge = {
        #     "debug": True,
        #     "genome": "27005/29/1",
        #     "output_name": "mergeGenome_temp",
        #     "annotations_to_merge": [
        #         {
        #             "annotation_source": ["KEGG KOs"],
        #             "annotation_weight": 0.9
        #         }, {
        #             "annotation_source": ["KOFAM93"],
        #             "annotation_weight": 0.5
        #         }, {
        #             "annotation_source": ['ProkkaAnnotation:2.1.5:EC:2020_10_13_20_07_38'],
        #             "annotation_weight": 1
        #         }
        #     ],
        #     "annotation_threshold": 1.0,
        #     "keep_best_annotation_only": 0,
        #     "description": "Merged annotations",
        #     "workspace_name": self.wsName
        # }
        # ret = self.serviceImpl.merge_metabolic_annotations(self.ctx, params_merge)
