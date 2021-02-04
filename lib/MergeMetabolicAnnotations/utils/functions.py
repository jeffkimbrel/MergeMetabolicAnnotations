import pandas as pd
import os
import logging
import yaml
import datetime
import json

import holoviews as hv
from holoviews import opts
from bokeh.models import HoverTool
hv.extension('bokeh')


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
    report.append("<style>* {font-family: sans-serif; font-size: 14px}</style>")

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


def get_event_lists(ontology):
    # print(type(ontology['feature_types']))

    gene_features = [k for k, v in ontology['feature_types'].items() if v == "gene"]

    events = {}
    for event in ontology["events"]:
        event_id = event['event_id']
        events[event_id] = {'genes': [],
                            'terms': [],
                            'msrxns': [],
                            'gene_msrxns': [],
                            'description': event['description'],
                            'timestamp': event['timestamp'],
                            'method': event['method'],
                            'method_version': event['method_version'],
                            'ontology_id': event['ontology_id']
                            }

        for gene in event["ontology_terms"]:
            if gene in gene_features:
                events[event_id]['genes'].append(gene)
                for entry in event["ontology_terms"][gene]:
                    if "term" in entry.keys():
                        events[event_id]['terms'].append(entry['term'])

                    if "modelseed_ids" in entry.keys():
                        events[event_id]['msrxns'] += entry['modelseed_ids']
                        for msrxn in entry['modelseed_ids']:
                            events[event_id]['gene_msrxns'].append(gene + '_' + msrxn)

        events[event_id]['genes'] = list(set(events[event_id]['genes']))
        events[event_id]['terms'] = list(set(events[event_id]['terms']))
        events[event_id]['msrxns'] = list(set(events[event_id]['msrxns']))
        events[event_id]['gene_msrxns'] = list(set(events[event_id]['gene_msrxns']))

    return events


def html_get_ontology_summary(event_summary, output_directory, to_highlight=[]):
    '''
    Only counts gene features, ignores cds
    The highlight part is a list of descriptions to have their rows highlighted.
    '''

    output_file = os.path.join(output_directory, "get_ontology_summary.html")

    report = html_header()
    report.append(f'<h3>Compare Annotations Summary</h3>')

    report.append(
        '<table cellspacing="0" cellpadding="3" border="1"><tr><th>Description</th><th>Timestamp</th><th>KBase App</th><th>Ontology</th><th>Genes</th><th>Unique Terms</th><th>Unique ModelSEED rxns</th><th>Unique Gene/ModelSEED rxn pairs</th></tr>')

    # get counts and add to new line of table
    for event in event_summary:
        if event_summary[event]["description"] in to_highlight:
            report.append(f'<tr style="background-color:#6CD075"><td>{event_summary[event]["description"].split(":")[0]}</td><td>{event_summary[event]["timestamp"]}</td><td>{event_summary[event]["method"]} v{event_summary[event]["method_version"]}</td><td>{event_summary[event]["ontology_id"]}</td><td>{len(event_summary[event]["genes"])}</td><td>{len(event_summary[event]["terms"])}</td><td>{len(event_summary[event]["msrxns"])}</td><td>{len(event_summary[event]["gene_msrxns"])}</td></tr>')
        else:
            report.append(f'<tr><td>{event_summary[event]["description"].split(":")[0]}</td><td>{event_summary[event]["timestamp"]}</td><td>{event_summary[event]["method"]} v{event_summary[event]["method_version"]}</td><td>{event_summary[event]["ontology_id"]}</td><td>{len(event_summary[event]["genes"])}</td><td>{len(event_summary[event]["terms"])}</td><td>{len(event_summary[event]["msrxns"])}</td><td>{len(event_summary[event]["gene_msrxns"])}</td></tr>')

    report.append('</table>')

    if len(to_highlight) > 0:
        report.append(
            f'<span style="background-color:#6CD075;font-size:12px"><i>* these highlighted rows were used for the merge</i></span>')

    print(output_file)

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
    elif workflow == "unique":
        selected_events = []

    for event in ontology["events"]:
        if event["description"] not in added_ontologies:  # keeps duplicates from being added twice
            if workflow == "unique":  # add all, don't filter
                ontology_selected['events'].append(event)
                added_ontologies.append(event["description"])
            else:
                if len(selected_events) == 0:  # if nothing is selected, then grab all events
                    if workflow == "merge":
                        event['annotation_weight'] = 1
                    ontology_selected['events'].append(event)
                    added_ontologies.append(event["description"])
                else:  # then grab only events in selected events, which is different for compare vs merge
                    if workflow == "compare":
                        if event["description"] in selected_events:
                            ontology_selected['events'].append(event)
                            added_ontologies.append(event["description"])
                    elif workflow == "merge":
                        for selected_event in selected_events:

                            # add if event in selected events, or if the first annotation_source is empty
                            if event["description"] in selected_event["annotation_source"] or len(selected_events[0]['annotation_source']) == 0:
                                event['annotation_weight'] = selected_event["annotation_weight"]
                                ontology_selected['events'].append(event)
                                added_ontologies.append(event["description"])

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

            # if best only is true and best_score is above threshold, use best_score as new threshold
            if best_score > params["annotation_threshold"]:
                gene_threshold = best_score
            else:
                gene_threshold = params["annotation_threshold"]
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


def plot_totals(event_summary, output_directory):

    df = pd.DataFrame(columns=['DESCRIPTION', 'GENES', 'TERMS', 'MSRXNS', 'GENE_MSRXNS'])

    for event in event_summary:
        df = df.append(pd.Series(data={
            'DESCRIPTION': event_summary[event]["description"],
            'GENES': len(event_summary[event]["genes"]),
            'TERMS': len(event_summary[event]["terms"]),
            'MSRXNS': len(event_summary[event]["msrxns"]),
            'GENE_MSRXNS': len(event_summary[event]["gene_msrxns"])}),
            ignore_index=True)

    df = pd.melt(df, id_vars=['DESCRIPTION'], var_name="TYPE", value_name="COUNT")

    def group_by(group, **kwargs):
        if group == "Description":
            bars = hv.Bars(df, kdims=['DESCRIPTION', 'TYPE'])
        elif group == "Type":
            bars = hv.Bars(df, kdims=['TYPE', 'DESCRIPTION'])

        return bars

    sets = ['Description', 'Type']
    bars = hv.DynamicMap(group_by, kdims='Group').redim.values(Group=sets)

    bars.opts(
        opts.Bars(color=hv.Cycle('Colorblind'), show_legend=False, stacked=False,
                  tools=['hover'], width=150*len(event_summary.keys()), height=600, xrotation=90))

    p_path = os.path.join(output_directory, 'totals.html')
    hv.output(widget_location='top')
    hv.save(bars, p_path, backend='bokeh')

    return {'path': output_directory,
            'name': os.path.basename(p_path),
            'description': 'Totals Report'}


def plot_agreements(event_summary, output_directory):

    def load_it_up(annotation_type, **kwargs):
        shared_genes = pd.DataFrame(columns=['A', 'B', 'AGREE', 'DISAGREE'])

        for event1 in event_summary:
            a = set(event_summary[event1][annotation_type])
            for event2 in event_summary:
                if event1 != event2:
                    b = set(event_summary[event2][annotation_type])
                    shared_genes = shared_genes.append(pd.Series(data={
                        'A': event1, 'B': event2, 'AGREE': len(a & b), 'DISAGREE': len(a - b)}), ignore_index=True)
                else:
                    shared_genes = shared_genes.append(pd.Series(data={
                        'A': event1, 'B': event2, 'AGREE': None, 'DISAGREE': None}), ignore_index=True)

        h1 = hv.HeatMap(shared_genes, ['A', 'B'], ['AGREE', 'DISAGREE'],
                        label="Agree (in row and column)")
        h2 = hv.HeatMap(shared_genes, ['A', 'B'], ['DISAGREE', 'AGREE'],
                        label="Disagree (in column only)")

        h1.opts(cmap='bgy')
        h2.opts(cmap='YlOrRd')

        return (h1 + h2).cols(1)

    sets = ['genes', 'msrxns', 'gene_msrxns']
    heatmap = hv.DynamicMap(load_it_up, kdims='Type').redim.values(Type=sets)

    hover = HoverTool(tooltips=[("Column", "@A"),
                                ("Row", "@B"),
                                ("Agree", "@AGREE"),
                                ("Disagree", '@DISAGREE')])

    heatmap.opts(
        opts.HeatMap(width=700,
                     height=100+75*len(event_summary.keys()),
                     tools=[hover],
                     axiswise=False,
                     logz=False,
                     invert_yaxis=True,
                     labelled=[],
                     toolbar='right',
                     colorbar=True,
                     xrotation=30),
        opts.VLine(line_color='black'))

    p_path = os.path.join(output_directory, 'agreements.html')
    hv.output(widget_location='top')
    hv.save(heatmap, p_path, backend='bokeh')

    return {'path': output_directory,
            'name': os.path.basename(p_path),
            'description': 'Agreements Report'}


def compare_report_stack(html_reports, event_summary, output_directory, to_highlight=[]):
    html_reports.append(html_get_ontology_summary(event_summary, output_directory, to_highlight))
    html_reports.append(plot_totals(event_summary, output_directory))
    html_reports.append(plot_agreements(event_summary, output_directory))
    html_reports.append(plot_csc(event_summary, output_directory))
    html_reports.append(plot_csc2(event_summary, output_directory))

    return html_reports


def plot_csc(event_summary, output_directory):

    def get_longest(type, modified_event_summary, aggregated_events):
        current_longest = {'description': '', 'a': 0}
        for event in modified_event_summary:
            a = len(list(set(modified_event_summary[event][type]) - set(aggregated_events)))
            if a >= current_longest['a']:  # need >= in case the last one is 0
                current_longest = {'description': event, 'a': a}
        return current_longest['description']

    def load_it_up(type, **kwargs):
        processed_events = {}
        aggregated_events = []
        original_event_count = len(event_summary.keys())
        modified_event_summary = event_summary.copy()
        bar_order = []

        df = pd.DataFrame(columns=['DESCRIPTION', 'COMPARISON', 'COUNT'])

        while len(processed_events.keys()) < original_event_count:
            l = modified_event_summary.pop(get_longest(
                type, modified_event_summary, aggregated_events))

            new_events = l[type]
            for event in processed_events:
                df = df.append(pd.Series(data={
                    'DESCRIPTION': l['description'], 'COMPARISON': event, 'COUNT': 0 - len(set(new_events) & set(processed_events[event]))}), ignore_index=True)
                new_events = set(new_events) - set(processed_events[event])

            processed_events[l['description']] = new_events
            aggregated_events = list(set(aggregated_events).union(set(new_events)))
            df = df.append(pd.Series(data={
                'DESCRIPTION': l['description'], 'COMPARISON': l['description'], 'COUNT': len(new_events)}), ignore_index=True)

        # create a line plot
        df_to_self = df[df['DESCRIPTION'] == df['COMPARISON']]
        df_to_self = df_to_self.drop(columns=['COMPARISON'])
        line_plot = hv.Curve(df_to_self).opts(color='k')
        dot_plot = hv.Scatter(df_to_self).opts(color='k', marker='+', size=10)

        bars = hv.Bars(df, ['DESCRIPTION', 'COMPARISON'], 'COUNT')
        bars = bars * line_plot * dot_plot

        return bars

    sets = ['genes', 'msrxns', 'gene_msrxns']
    bars = hv.DynamicMap(load_it_up, kdims='Type').redim.values(Type=sets)

    bars.opts(
        opts.Bars(color=hv.Cycle('Colorblind'), invert_axes=True, show_legend=False, stacked=True,
                  tools=['hover'], width=1000, height=100+50*len(event_summary.keys()), xrotation=90))

    p_path = os.path.join(output_directory, 'csc.html')
    hv.output(widget_location='top')
    hv.save(bars, p_path, backend='bokeh')

    return {'path': output_directory,
            'name': os.path.basename(p_path),
            'description': 'CSC Report'}


def plot_csc2(event_summary, output_directory, descript_truncate=50):

    def get_longest(type, modified_event_summary, aggregated_events):
        current_longest = {'description': '', 'a': 0}
        for event in modified_event_summary:
            a = len(list(set(modified_event_summary[event][type]) - set(aggregated_events)))
            if a >= current_longest['a']:  # need >= in case the last one is 0
                current_longest = {'description': event, 'a': a}
        return current_longest['description']

    def load_it_up(type, **kwargs):
        processed_events = {}
        aggregated_events = []
        original_event_count = len(event_summary.keys())
        modified_event_summary = event_summary.copy()
        # bar_order = []
        baseline = 0

        df = pd.DataFrame(columns=['DESCRIPTION', 'COMPARISON', 'COUNT'])

        while len(processed_events.keys()) < original_event_count:
            l = modified_event_summary.pop(get_longest(
                type, modified_event_summary, aggregated_events))

            new_events = l[type]
            sub_baseline = baseline

            # get what overlaps with already processed events
            for event in processed_events:
                overlap = set(new_events) & set(processed_events[event])
                new_events = set(new_events) - set(processed_events[event])
                df = df.append(pd.Series(data={
                    'DESCRIPTION': l['description'][:descript_truncate], 'COMPARISON': event[:descript_truncate], 'COUNT': len(overlap), 'HIGH': sub_baseline, 'LOW': sub_baseline-len(overlap)}), ignore_index=True)
                sub_baseline -= len(overlap)

            processed_events[l['description']] = new_events
            aggregated_events = list(set(aggregated_events).union(set(new_events)))

            # if anything is left, add it has new events
            df = df.append(pd.Series(data={
                'DESCRIPTION': l['description'][:descript_truncate], 'COMPARISON': l['description'][:descript_truncate], 'COUNT': len(new_events), 'HIGH': baseline + len(new_events), 'LOW': baseline}), ignore_index=True)
            baseline = len(aggregated_events)

        print(df)

        # # flip the df order so plot order is flipped.
        # df = df.loc[::-1].reset_index(drop=True)
        # print(df)

        seg = hv.Segments(df, [hv.Dimension('LOW', label='Count'),
                               hv.Dimension('DESCRIPTION', label='Genome Event'),
                               'HIGH',
                               'DESCRIPTION'])

        return seg

    sets = ['genes', 'msrxns', 'gene_msrxns']
    seg = hv.DynamicMap(load_it_up, kdims='Type').redim.values(Type=sets)

    csc_hover = HoverTool(tooltips=[("Main", "@DESCRIPTION"),
                                    ("Comp", "@COMPARISON"),
                                    ("Count", "@COUNT")])

    seg.opts(line_width=40, color='COMPARISON', cmap='bgy', height=100+50*len(event_summary.keys()),
             width=800, tools=[csc_hover], invert_axes=False, xrotation=90)

    # bars.opts(
    #     opts.Bars(color=hv.Cycle('Colorblind'), invert_axes=True, show_legend=False, stacked=True,
    #               tools=['hover'], width=1000, height=100+50*len(event_summary.keys()), xrotation=90))

    p_path = os.path.join(output_directory, 'csc2.html')
    hv.output(widget_location='top')
    hv.save(seg, p_path, backend='bokeh')

    return {'path': output_directory,
            'name': os.path.basename(p_path),
            'description': 'CSC Report'}


if __name__ == "__main__":
    ontology_selected = json.loads(
        open("/Users/kimbrel1/Desktop/get_ontology_dump_after_merge.json", "r").read())

    d = get_event_lists(ontology_selected)
    p_path = plot_csc2(d, "/Users/kimbrel1/Desktop/")
    print(p_path)
