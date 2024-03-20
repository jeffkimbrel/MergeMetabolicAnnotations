#-*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os

from .utils.ImportAnnotationsUtil import ImportAnnotationsUtil
from MergeMetabolicAnnotations.utils.ImportAnnotationsUtil import ImportAnnotationsUtil

from .utils.ImportBulkAnnotationsUtil import ImportBulkAnnotationsUtil
from MergeMetabolicAnnotations.utils.ImportBulkAnnotationsUtil import ImportBulkAnnotationsUtil

from .utils.CompareAnnotationsUtil import CompareAnnotationsUtil
from MergeMetabolicAnnotations.utils.CompareAnnotationsUtil import CompareAnnotationsUtil

from .utils.MergeAnnotationsUtil import MergeAnnotationsUtil
from MergeMetabolicAnnotations.utils.MergeAnnotationsUtil import MergeAnnotationsUtil

from installed_clients.KBaseReportClient import KBaseReport
#END_HEADER


class MergeMetabolicAnnotations:
    '''
    Module Name:
    MergeMetabolicAnnotations

    Module Description:
    A KBase module: MergeMetabolicAnnotations
    This module implements tools for importing, comparing and merging
    3rd party metabolic annotations.
    '''

    ########WARNING FOR GEVENT USERS #######noqa
    #Since asynchronous IO can lead to methods - even the same method -
    #interrupting each other, you must be *very* careful when using global
    #state. A method could easily clobber the state set by another while
    #the latter method is running.
    #########################################noqa
    VERSION = "0.0.1"
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    #config contains contents of config file in a hash or None if it couldn't
    #be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        self.config = config
        self.config['SDK_CALLBACK_URL'] = self.callback_url
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass

    def import_annotations(self, ctx, params):
        """
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> unspecified object
        """
        #ctx is the context object
        #return variables are: output
        #BEGIN import_annotations
        import_runner = ImportAnnotationsUtil(self.config)
        output = import_runner.run(ctx, params)
        #END import_annotations

        print(output)

        #At some point might do deeper type checking...
        if not isinstance(output, object):
            raise ValueError('Method import_annotations return value ' +
                             'output is not type object as required.')
        #return the results
        return [output]

    def import_bulk_annotations(self, ctx, params):
        """
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> unspecified object
        """
        #ctx is the context object
        #return variables are: output
        #BEGIN import_bulk_annotations
        import_runner = ImportBulkAnnotationsUtil(self.config)
        output = import_runner.run(ctx, params)
        #END import_bulk_annotations

        #At some point might do deeper type checking...
        if not isinstance(output, object):
            raise ValueError('Method import_bulk_annotations return value ' +
                             'output is not type object as required.')
        #return the results
        return [output]

    def compare_metabolic_annotations(self, ctx, params):
        """
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> unspecified object
        """
        #ctx is the context object
        #return variables are: output
        #BEGIN compare_metabolic_annotations
        import_runner = CompareAnnotationsUtil(self.config)
        output = import_runner.run(ctx, params)
        #END compare_metabolic_annotations

        #At some point might do deeper type checking...
        if not isinstance(output, object):
            raise ValueError('Method compare_metabolic_annotations return value ' +
                             'output is not type object as required.')
        #return the results
        return [output]

    def merge_metabolic_annotations(self, ctx, params):
        """
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> unspecified object
        """
        #ctx is the context object
        #return variables are: output
        #BEGIN merge_metabolic_annotations
        import_runner = MergeAnnotationsUtil(self.config)
        output = import_runner.run(ctx, params)
        #END merge_metabolic_annotations

        #At some point might do deeper type checking...
        if not isinstance(output, object):
            raise ValueError('Method merge_metabolic_annotations return value ' +
                             'output is not type object as required.')
        #return the results
        return [output]

    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
