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
        # params_import = {
        #     "debug": True,
        #     "ontology": "KO",
        #     "annotation_file": "20201123_blastkoala_forUpload.txt",
        #     "description": "import_test",
        #     "genome": "52279/9/1",
        #     "output_name": "import_genome",
        #     "workspace_name": self.wsName
        # }
        # ret = self.serviceImpl.import_annotations(self.ctx, params_import)

        # bulk import app
        # params_import = {
        #     "debug": True,
        #     "annotation_file": "pf5_bulk_test2.txt",
        #     "genome": "52279/9/1",
        #     "output_name": "import_bulk_genome",
        #     "workspace_name": self.wsName
        # }
        # ret = self.serviceImpl.import_bulk_annotations(self.ctx, params_import)

        # compare app
        params_compare = {
            "debug": True,
            "genome": "116449/93/2",
            "output_name": "compare_test",
            "annotations_to_compare": [],
            # "annotations_to_compare": ["ProkkaAnnotation:3.2.1:EC:2021_01_14_16_28_25", "koalaFams_v96_KO:2021_01_14_16_32_51"],
            "workspace_name": self.wsName
        }
        ret = self.serviceImpl.compare_metabolic_annotations(self.ctx, params_compare)

        # merge app
        # params_merge = {
        #     "debug": True,
        #     "genome": "52279/12/1",
        #     "output_name": "merge_test",
        #     # "annotations_to_merge": [{
        #     #     "annotation_source": ["ProkkaAnnotation:3.2.1:EC:2021_01_14_16_28_25"],
        #     #     "annotation_weight": 0.6
        #     # }, {
        #     #     "annotation_source": ["BlastKOALA_V96_KOs:2021_01_21_02_23_59"],
        #     #     "annotation_weight": 0.6
        #     # }],
        #     "annotations_to_merge": [{
        #         "annotation_source": [],
        #         "annotation_weight": 1
        #     }],
        #     "annotation_threshold": 2,
        #     "keep_best_annotation_only": 0,
        #     "description": "merge t=1 best=F",
        #     "workspace_name": self.wsName
        # }
        # ret = self.serviceImpl.merge_metabolic_annotations(self.ctx, params_merge)
