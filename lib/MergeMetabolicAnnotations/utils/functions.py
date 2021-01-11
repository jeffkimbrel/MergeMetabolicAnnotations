import pandas as pd
import os
import logging
import yaml
import datetime


def df_to_ontology(params, staging_dir, pass_df=None):
    '''
    Takes the text file from staging, or the pandas df passed from the merge
    app, and converts to an ontology dictionary suitable from the annotation
    ontology API add_annotation_ontology_events() method
    '''

    if isinstance(pass_df, pd.DataFrame):
        annotations = pass_df
        method = "Merge Annotations"
    else:
        if 'debug' in params and params['debug'] is True:
            annotations_file_path = os.path.join(
                '/kb/module/test/test_data', params['annotation_file'])
        else:
            annotations_file_path = os.path.join(staging_dir, params['annotation_file'])

        annotations = pd.read_csv(annotations_file_path,
                                  sep='\t',
                                  header=None,
                                  names=['gene', 'term']
                                  )

        method = "Import Annotations"

    # remove duplicate rows, if any
    annotations = annotations.drop_duplicates()

    ontology = {
        'event_id': params['description'],
        'ontology_id': params['ontology'],
        'method': method,  # from above
        'method_version': get_app_version(),
        "timestamp": datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S"),
        'ontology_terms': {}
    }

    # add imported terms
    for index, row in annotations.iterrows():
        if pd.notnull(row['term']):
            if row['gene'] in ontology['ontology_terms']:
                ontology['ontology_terms'][row['gene']].append(
                    {'term': row['term']}
                )
            else:
                ontology['ontology_terms'][row['gene']] = [
                    {'term': row['term']}
                ]

    return ontology


def get_app_version():
    with open("/kb/module/kbase.yml", 'r') as stream:
        data_loaded = yaml.load(stream)
    return str(data_loaded['module-version'])
