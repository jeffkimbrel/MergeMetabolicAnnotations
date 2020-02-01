/*
A KBase module: MergeMetabolicAnnotations
This module implements tools for importing, comparing and merging 3rd party metabolic annotations.
*/

module MergeMetabolicAnnotations {

  typedef UnspecifiedObject ReportResults;

  funcdef import_annotations(mapping<string,UnspecifiedObject> params)
    returns (ReportResults output) authentication required;

  funcdef import_bulk_annotations(mapping<string,UnspecifiedObject> params)
    returns (ReportResults output) authentication required;

  funcdef compare_metabolic_annotations(mapping<string,UnspecifiedObject> params)
    returns (ReportResults output) authentication required;

  funcdef merge_metabolic_annotations(mapping<string,UnspecifiedObject> params)
    returns (ReportResults output) authentication required;

  funcdef run_MergeMetabolicAnnotations(mapping<string,UnspecifiedObject> params)
    returns (ReportResults output) authentication required;

};
