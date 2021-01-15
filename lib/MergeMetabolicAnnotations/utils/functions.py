import pandas as pd
import os
import logging
import yaml
import datetime


def df_to_ontology(params, pass_df=None):
    '''
    Takes the text file from staging, or the pandas df passed from the merge
    app, and converts to an ontology dictionary suitable from the annotation
    ontology API add_annotation_ontology_events() method

    The merge app also calls this, and it can use the same params that the
    import gives... all shared features are in both (except for the
    annotations_file which the html_add_ontology_summary needs to fix)
    '''

    if isinstance(pass_df, pd.DataFrame):
        annotations = pass_df
        method = "Merge Annotations"
    else:
        if 'debug' in params and params['debug'] is True:
            annotations_file_path = os.path.join(
                '/kb/module/test/test_data', params['annotation_file'])
        else:
            annotations_file_path = os.path.join("/staging/", params['annotation_file'])

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
        'description': params['description'],
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


def html_header():
    report = []
    report.append("<style>* {font-family: sans-serif;}</style>")

    return report


def html_add_ontology_summary(params, ontology, api_results, output_directory):

    output_file = os.path.join(output_directory, "add_ontology_summary.html")

    # Make report directory and copy over files
    report = html_header()

    report.append(f'<h3>Import Annotations Summary</h3>')

    report.append(f'<b>Import App version:</b> {ontology["method_version"]}<br>')
    report.append(f'<b>Timestamp:</b> {ontology["timestamp"]}<br>')

    # running this code with the merge app doesn't have an annotations file
    if "annotation_file" in params:
        report.append(f'<b>Annotations file:</b> {params["annotation_file"]}<br>')

    report.append(f'<b>Ontology ID:</b> {params["ontology"]}<br>')
    report.append(f'<b>Description:</b> {params["description"]}<br>')

    report.append(f'<b>Input Ref:</b> {params["genome"]}<br>')
    report.append(f'<b>Output Ref:</b> {api_results["output_ref"]}<br><br>')

    report.append(f'<b>Features in annotations file:</b> {len(ontology["ontology_terms"])}<br>')
    report.append(f'<b>Features (found):</b> {api_results["ftrs_found"]}<br>')
    report.append(f'<b>Features (not found):</b> {len(api_results["ftrs_not_found"])}<br>')

    if len(api_results["ftrs_not_found"]) > 0:
        report.append(
            f'<br><b>These genes were not found in the genome:</b> <br>{("<br>").join(api_results["ftrs_not_found"])}<br>')

    # Write to file
    with open(output_file, 'w') as f:
        for line in report:
            f.write(line + "\n")

    return {'path': output_directory,
            'name': os.path.basename(output_file),
            'description': 'HTML report for import_annotations app'}


def html_get_ontology_summary(ontology, output_directory):
    '''
    Only counts gene features, ignores cds
    '''

    output_file = os.path.join(output_directory, "get_ontology_summary.html")

    report = html_header()
    report.append(f'<h3>Compare Annotations Summary</h3>')

    report.append(
        '<table cellspacing="0" cellpadding="3" border="1"><tr><th>Description</th><th>Timestamp</th><th>App Version</th><th>Ontology</th><th>Genes</th><th>Unique Terms</th><th>Unique ModelSEED rxns</th></tr>')

    # get counts and add to new line of table
    gene_features = [k for k, v in ontology['feature_types'].items() if v == "gene"]

    for event in ontology["events"]:
        gene_count = 0

        terms = []
        msrxns = []

        for gene in event["ontology_terms"]:
            if gene in gene_features:
                gene_count += 1
                for entry in event["ontology_terms"][gene]:
                    if "term" in entry.keys():
                        terms.append(entry['term'])

                    if "modelseed_ids" in entry.keys():
                        msrxns += entry['modelseed_ids']

        # term_count = 0
        # msrxn_count = 0
        term_count = len(list(set(terms)))
        msrxn_count = len(list(set(msrxns)))

        report.append(
            f'<tr><td>{event["description"].split(":")[0]}</td><td>{event["timestamp"]}</td><td>{event["method_version"]}</td><td>{event["ontology_id"]}</td><td>{gene_count}</td><td>{term_count}</td><td>{msrxn_count}</td></tr>')

    report.append('</table>')

    # Write to file
    with open(output_file, 'w') as f:
        for line in report:
            f.write(line + "\n")

    return {'path': output_directory,
            'name': os.path.basename(output_file),
            'description': 'Summary Report'}


def filter_selected_ontologies(ontology, params, workflow="compare"):
    '''
    unique will not use params and just give all unique event_ids.
    compare and merge workflows filter out the ontologies to those selected in
    the UI, but the merge workflow also adds the annotation_weights to the
    events. both return all unique if no events are selected, and defaulting to
    a weight of 1 for the merge

    The unique functionality is added because of a current bug
    '''

    ontology_selected = {"events": [],
                         "feature_types": ontology["feature_types"]}

    added_ontologies = []

    # the list of selections have different names depending on the app
    if workflow == "compare":
        selected_events = params['annotations_to_compare']
    elif workflow == "merge":
        selected_events = params['annotations_to_merge']

    for event in ontology["events"]:
        if event["event_id"] not in added_ontologies:
            if workflow == "unique":
                ontology_selected['events'].append(event)
                added_ontologies.append(event["event_id"])
            else:
                if len(selected_events) == 0:
                    if workflow == "merge":
                        event['annotation_weight'] = 1
                    ontology_selected['events'].append(event)
                    added_ontologies.append(event["event_id"])
                else:
                    for selected_event in selected_events:
                        if event["event_id"] == selected_event["annotation_source"][0]:
                            if workflow == "merge":
                                event['annotation_weight'] = selected_event["annotation_weight"]
                            ontology_selected['events'].append(event)
                            added_ontologies.append(event["event_id"])

    return ontology_selected


def merge_ontology_events(ontology):
    '''
    The annotation ontology api can put annotations in the cds features as well
    as the gene features. This code only considers gene features, and ignores
    annotations in cds features. I think they only get added to cds features if
    they were aliases
    '''

    # get counts and add to new line of table
    gene_features = [k for k, v in ontology['feature_types'].items() if v == "gene"]

    ontology_merged = {}

    for event in ontology["events"]:
        for gene in event["ontology_terms"]:
            if gene in gene_features:
                for entry in event["ontology_terms"][gene]:
                    if "modelseed_ids" in entry.keys():
                        for MSRXN in entry['modelseed_ids']:
                            if gene in ontology_merged:
                                if MSRXN in ontology_merged[gene]:
                                    ontology_merged[gene][MSRXN] += event['annotation_weight']
                                else:
                                    ontology_merged[gene][MSRXN] = event['annotation_weight']
                            else:
                                ontology_merged[gene] = {MSRXN: event['annotation_weight']}

    return ontology_merged


def score_mergers(ontology_merged, params):
    '''
    returns a pandas dataframe suitable for the import annotations workflow
    '''

    df = pd.DataFrame(columns=['gene', 'term', 'score', 'pass'])

    for gene_id in ontology_merged:
        if params["keep_best_annotation_only"] == 1:
            best_score = 0
            for MSRXN in ontology_merged[gene_id]:
                if ontology_merged[gene_id][MSRXN] > best_score:
                    best_score = ontology_merged[gene_id][MSRXN]

            # replace threshold value with best score, if "best only" is true
            # righ tnow, this happens whether or not the best_score is above the
            # threshold, but may need to only happen if best_score > treshold
            gene_threshold = best_score

        else:
            gene_threshold = params["annotation_threshold"]

        for MSRXN in ontology_merged[gene_id]:
            if ontology_merged[gene_id][MSRXN] >= gene_threshold:
                df = df.append(pd.Series(data={
                               'gene': gene_id, 'term': MSRXN, 'score': ontology_merged[gene_id][MSRXN], 'pass': 1}), ignore_index=True)

            else:
                df = df.append(pd.Series(data={
                               'gene': gene_id, 'term': MSRXN, 'score': ontology_merged[gene_id][MSRXN], 'pass': 0}), ignore_index=True)

    # returns all results
    return df
