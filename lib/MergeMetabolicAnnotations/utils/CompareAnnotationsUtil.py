import os
import datetime
import logging
import json
import uuid
import re
import pandas as pd

from bokeh.plotting import figure, output_file
from bokeh.io import save
from bokeh.palettes import inferno, viridis
from bokeh.models import ColumnDataSource, FactorRange
from bokeh.transform import factor_cmap
from bokeh.models import HoverTool

from installed_clients.GenomeAnnotationAPIClient import GenomeAnnotationAPI
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.WorkspaceClient import Workspace as Workspace
from installed_clients.KBaseReportClient import KBaseReport

import MergeMetabolicAnnotations.utils.utils as mu


class CompareAnnotationsUtil:

    workdir = 'tmp/work/'
    staging_dir = "/staging/"
    datadir = "/kb/module/data/"

    # to make case-insensitive, these are upper case, and the utils.get_ontology_dict() will convert to upper
    translation_locations = {'EC': 'EBI_EC.ModelSEED.json',
                             'RO': 'KEGG_RXN.ModelSEED.json',
                             'KO': 'KEGG_KO.ModelSEED.json',
                             'META': 'Metacyc_RXN.ModelSEED.json',
                             'SSO': 'SSO.ModelSEED.json',
                             'GO': 'GO.ModelSEED.json'}

    legacy_codes = {'MODELSEED': 'MSRXN',
                    'KEGGKO': 'KO',
                    'KEGGRO': 'RO',
                    'METACYC': 'META'
                    }

    def __init__(self, config):
        os.makedirs(self.workdir, exist_ok=True)
        self.config = config
        self.timestamp = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.genome_api = GenomeAnnotationAPI(self.callback_url)
        self.dfu = DataFileUtil(self.callback_url)
        self.gfu = GenomeFileUtil(self.callback_url)
        self.kbr = KBaseReport(self.callback_url)
        self.ws_client = Workspace(config["workspace-url"])

        self.events = {}

    def get_ontology_events(self, params):

        if 'ontology_events' in self.genome:
            for event, ontology in enumerate(self.genome['ontology_events']):
                if 'description' not in ontology:
                    ontology['description'] = ontology['method']
                if ontology['description'] in params['annotations_to_compare'] or len(params['annotations_to_compare']) == 0:
                    self.events[event] = {}

                    ontology["id"] = ontology["id"].upper()

                    if ontology["id"] in self.legacy_codes:
                        ontology["id"] = self.legacy_codes[ontology["id"]]

                    for term in ontology:
                        self.events[event][term] = ontology[term]
        else:
            logging.info("No ontology events in this genome!")

        logging.info(self.events)

    def get_translations(self):
        self.translations = {}

        for ontology_type in self.translation_locations:
            translations_path = self.datadir + "/" + self.translation_locations[ontology_type]
            ontology_translations = json.loads(open(translations_path, "r").read())
            self.translations[ontology_type] = {}

            for term in ontology_translations['translation']:

                for entry in ontology_translations['translation'][term]['equiv_terms']:
                    if entry['equiv_term'] != None:
                        rxn = entry['equiv_term']
                        # logging.info(rxn)
                        if not rxn.upper().startswith('MSRXN:'):
                            rxn = 'MSRXN:' + rxn
                        # logging.info(rxn)

                        # Add ontology type if not present
                        if not term.upper().startswith(ontology_type + ':'):
                            term = ontology_type + ":" + term

                        self.translations[ontology_type][term] = rxn

    def search_for_ec(self, line):
        ecList = re.findall(r"\(*[0-9]+\.[0-9\-]+\.[0-9\-]+\.[0-9\-]+", line)
        return(ecList)

    def summarize_gto(self, params):
        summary = {"genes": {},
                   "terms": {},
                   "rxns": {},
                   "ontology_events": {},
                   "orphan_terms": {}
                   }

        # add ontology events
        for ontology_event in self.events:
            summary['ontology_events'][ontology_event] = self.events[ontology_event]

        # add gene id to summary
        for feature in self.genome['features']:
            gene_id = feature['id']
            summary["genes"][gene_id] = {"terms": {},
                                         "rxns": {}
                                         }

            # get ontology term
            if "ontology_terms" in feature:

                for ontology_term_type in feature['ontology_terms']:
                    '''
                    note, ontology_term_type might be a legacy term, and will need to be converted later after making the term_dict
                    '''

                    # ontology_term_type = ontology_term_type.upper()
                    # logging.info(ontology_term_type)
                    # if ontology_term_type in self.legacy_codes:
                    #     ontology_term_type = self.legacy_codes[ontology_term_type]
                    # logging.info(ontology_term_type)

                    term_dict = feature['ontology_terms'][ontology_term_type]

                    for term in term_dict:
                        for ontology_event in term_dict[term]:

                            # is this ontology event in the user-selected list?
                            if ontology_event in self.events:

                                rxn = "none"

                                # get rxn, convert to upper case to make case insensitive
                                ontology_type = summary['ontology_events'][ontology_event]['id'].upper(
                                )

                                # fix annotation term to fit with style
                                if not term.startswith(ontology_type):
                                    if term.startswith('ec:'):
                                        term = term.upper()
                                    else:
                                        term = ontology_type + ':' + term

                                # convert terms to rxns
                                if term in self.translations[ontology_type]:
                                    rxn = self.translations[ontology_type][term]
                                else:
                                    if ontology_event in summary["orphan_terms"]:
                                        summary["orphan_terms"][ontology_event].append(term)
                                        summary["orphan_terms"][ontology_event] = list(
                                            set(summary["orphan_terms"][ontology_event]))
                                    else:
                                        summary["orphan_terms"][ontology_event] = [term]

                                # terms
                                if term in summary["genes"][gene_id]['terms']:
                                    summary["genes"][gene_id]['terms'][term].append(ontology_event)
                                else:
                                    summary["genes"][gene_id]['terms'][term] = [ontology_event]

                                if term in summary['terms']:
                                    summary['terms'][term].append(ontology_event)
                                    summary['terms'][term] = list(set(summary['terms'][term]))
                                else:
                                    summary['terms'][term] = [ontology_event]

                                # rxns
                                if rxn != "none":
                                    if rxn in summary["genes"][gene_id]['rxns']:
                                        summary["genes"][gene_id]['rxns'][rxn].append(
                                            ontology_event)
                                    else:
                                        summary["genes"][gene_id]['rxns'][rxn] = [ontology_event]

                                    if rxn in summary['rxns']:
                                        summary['rxns'][rxn].append(ontology_event)
                                        summary['rxns'][rxn] = list(set(summary['rxns'][rxn]))
                                    else:
                                        summary['rxns'][rxn] = [ontology_event]

        with open(os.path.join(self.scratch, "summary_dump.json"), 'w') as outfile:
            json.dump(summary, outfile, indent=2)

        return summary

    def html_summary(self, params, summary):

        # convert gto summary for this report
        html_summary_report = {}

        for ontology_event in summary['ontology_events']:
            html_summary_report[ontology_event] = {"gene": [], "term": [], "rxn": []}

        for gene in summary["genes"]:
            for term in summary["genes"][gene]['terms']:
                for ontology_event in summary["genes"][gene]['terms'][term]:
                    html_summary_report[ontology_event]['gene'].append(gene)
                    html_summary_report[ontology_event]['term'].append(term)

                    html_summary_report[ontology_event]['gene'] = list(
                        set(html_summary_report[ontology_event]['gene']))
                    html_summary_report[ontology_event]['term'] = list(
                        set(html_summary_report[ontology_event]['term']))

            for rxn in summary["genes"][gene]['rxns']:
                for ontology_event in summary["genes"][gene]['rxns'][rxn]:
                    html_summary_report[ontology_event]['rxn'].append(rxn)
                    html_summary_report[ontology_event]['gene'].append(gene)

                    html_summary_report[ontology_event]['rxn'] = list(
                        set(html_summary_report[ontology_event]['rxn']))
                    html_summary_report[ontology_event]['gene'] = list(
                        set(html_summary_report[ontology_event]['gene']))

        output_html_files = list()

        # Make report directory and copy over files
        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        os.mkdir(output_directory)
        result_file_path = os.path.join(output_directory, 'compare_annotations_summary.html')

        # make html
        table_lines = []
        table_lines.append(f'<h2>Compare Annotations</h2>')
        table_lines.append(f'<h3>Summary</h3>')
        table_lines.append(
            '<table cellspacing="0" cellpadding="3" border="1"><tr><th>EVENT</th><th>DESCRIPTION</th><th>TYPE</th><th>GENES</th><th>TERMS</th><th>RXNS</th></tr>')
        for event in sorted(html_summary_report.keys()):
            # RAST/PROKKA don't have descriptions, but they have methods
            description = self.events[event].get('description', self.events[event]['method'])
            type = self.events[event]['id']
            genes_list = html_summary_report[event]['gene']
            terms_list = html_summary_report[event]['term']
            rxns_list = html_summary_report[event]['rxn']
            table_lines.append('<tr><td>' + str(event) + '</td><td>' + description + '</td><td>' + type + '</td><td>' + str(
                len(set(genes_list))) + '</td><td>' + str(len(terms_list)) + '</td><td>' + str(len(rxns_list)) + '</td></tr>')
        table_lines.append('</table>')

        # Write to file
        with open(result_file_path, 'w') as result_file:
            for line in table_lines:
                result_file.write(line + "\n")

        output_html_files.append(
            {'path': output_directory,
             'name': os.path.basename(result_file_path),
             'description': 'Summary Report'})

        # bokeh plots
        totals_file_path = os.path.join(output_directory, 'totals.html')
        output_file(totals_file_path, title="Totals")
        totals = self.plot_totals(summary)
        save(totals)
        output_html_files.append(
            {'path': output_directory,
             'name': os.path.basename(totals_file_path),
             'description': 'Ontology Totals'})

        csc_file_path = os.path.join(output_directory, 'csc.html')
        output_file(csc_file_path, title="CSC")
        csc = self.plot_csc2(summary)
        save(csc)
        output_html_files.append(
            {'path': output_directory,
             'name': os.path.basename(csc_file_path),
             'description': 'Cumulative Sum Plot'})

        # finalize html reports
        report_params = {
            'message': '',
            'html_links': output_html_files,
            'direct_html_link_index': 0,
            'workspace_name': params['workspace_name'],
            'report_object_name': f'compare_annotations_{uuid.uuid4()}'}

        output = self.kbr.create_extended_report(report_params)

        return {'report_name': output['name'],
                'report_ref': output['ref']}

# plotting functions
    def plot_totals(self, summary):
        descriptions = {}
        for o in summary["ontology_events"]:
            descriptions[o] = summary["ontology_events"][o].get(
                'description', summary["ontology_events"][o]['method']) + '_' + str(o)
            logging.info(descriptions[o])

        totals = {}
        for event in summary['ontology_events'].keys():
            totals[str(event)] = {'genes': [],
                                  'rxns': [],
                                  'terms': []}

        # genes
        for gene in summary['genes']:
            for term in summary['genes'][gene]['terms']:
                for event in summary['genes'][gene]['terms'][term]:
                    totals[str(event)]['genes'].append(gene)

        # terms
        for term in summary['terms']:
            for event in summary['terms'][term]:
                totals[str(event)]['terms'].append(term)

        # rxns
        for rxn in summary['rxns']:
            for event in summary['rxns'][rxn]:
                totals[str(event)]['rxns'].append(rxn)

        # sums
        events = []
        types = ['genes', 'terms', 'rxns']

        gene_counts = []
        rxn_counts = []
        term_counts = []

        for event in totals:
            logging.info(event)
            events.append(descriptions[int(event)])
            gene_counts.append(len(set(totals[event]['genes'])))
            rxn_counts.append(len(set(totals[event]['rxns'])))
            term_counts.append(len(set(totals[event]['terms'])))

        data = {'events': events,
                'genes': gene_counts,
                'terms': term_counts,
                'rxns': rxn_counts
                }

        x = [(event, type) for event in events for type in types]

        counts = sum(zip(data['genes'], data['terms'], data['rxns']), ())
        source = ColumnDataSource(data=dict(x=x, counts=counts))

        p = figure(y_range=FactorRange(*x),
                   plot_height=400,
                   plot_width=1000,
                   title="Unique Counts per Annotation Event",
                   tools="wheel_zoom,box_zoom,reset,save")

        p.hbar(y='x',
               right='counts',
               height=0.9,
               source=source,
               line_color="black",
               fill_color=factor_cmap('x',
                                      palette=inferno(len(types)),
                                      factors=types,
                                      start=1,
                                      end=2))

        p.x_range.start = 0
        p.y_range.range_padding = 0.1
        p.yaxis.major_label_orientation = "horizontal"
        p.yaxis.subgroup_label_orientation = "horizontal"
        p.yaxis.group_label_orientation = "horizontal"
        p.ygrid.grid_line_color = None
        p.title.text_font_size = '12pt'
        p.xaxis.major_label_text_font_size = "12pt"
        p.yaxis.major_label_text_font_size = "12pt"
        p.yaxis.group_text_font_size = "12pt"
        p.add_tools(HoverTool(tooltips=[("Type", "@x"), ("Count", "@counts")]))

        return(p)

    def plot_csc2(self, summary, summary_type="rxns"):
        descriptions = {}
        for o in summary["ontology_events"]:
            descriptions[o] = summary["ontology_events"][o].get(
                'description', summary["ontology_events"][o]['method']) + ' (' + summary["ontology_events"][o]['id'] + ' #' + str(o) + ')'

        events = sorted(summary['ontology_events'].keys())
        rxns = summary[summary_type]

        # convert to sets
        rxns_in_events = dict((int(el), set()) for el in events)
        for rxn in rxns:
            for event in rxns[rxn]:
                rxns_in_events[event].add(rxn)

        winning_sets = {}
        winning_order = []
        baseline = 0
        df = pd.DataFrame(columns=["E", "C", "T", "L", "R"])
        # E=event, C=comparison, T=total, L=left, R=right

        for _ in range(len(rxns_in_events)):

            current_right = baseline
            current_left = baseline

            # get current winner
            longest_set_key = self.longest_set(rxns_in_events, winning_sets)

            # compare current winner to all past winners
            current = rxns_in_events[longest_set_key]
            for past_winner in winning_order:
                overlap = len(winning_sets[past_winner] & current)
                current_left -= overlap
                row = [descriptions[longest_set_key],  # E
                       descriptions[past_winner],  # C
                       overlap,  # T
                       current_left,  # L
                       current_left + overlap]  # R
                df.loc[len(df)] = row
                current = current - winning_sets[past_winner]

            # process current winner
            row = [descriptions[longest_set_key],  # E
                   descriptions[longest_set_key],  # C
                   len(current),  # T
                   current_right,  # L
                   current_right + len(current)]  # R

            df.loc[len(df)] = row  # add to df
            baseline += len(current)

            # move current winner to past winners
            winning_sets[longest_set_key] = rxns_in_events[longest_set_key]
            winning_order.append(longest_set_key)
            rxns_in_events[longest_set_key] = set()

        source = ColumnDataSource(df)

        type1_colormap = factor_cmap('E', palette=viridis(
            len(df.E.unique())), factors=df.E.unique())
        type2_colormap = factor_cmap('C', palette=viridis(
            len(df.C.unique())), factors=df.C.unique())

        p = figure(y_range=df.E.unique().tolist()[::-1],  # .tolist()[::-1] reverses the list.
                   plot_height=300,
                   plot_width=1000,
                   title="Annotation events ranked by \'" + str(summary_type) + "\' contribution",
                   tools="wheel_zoom,box_zoom,reset,save")

        p.hbar(y='E',
               height=0.9,
               left='L',
               right='R',
               source=source,
               fill_color=type2_colormap,
               line_color="black")

        p.add_tools(HoverTool(tooltips=[("Total", "@T"), ("Comparison", "@C")]))
        p.title.text_font_size = '12pt'
        p.xaxis.major_label_text_font_size = "12pt"
        p.yaxis.major_label_text_font_size = "12pt"
        return p

    def longest_set(self, s, w):
        s = s.copy()
        for event in s:
            for winner in w:
                s[event] = s[event] - w[winner]

        # https://stackoverflow.com/a/21839239
        max_key, max_value = max(s.items(), key=lambda x: len(x[1]))
        return(max_key)

    def run(self, ctx, params):

        # collect some metadata
        self.genome = mu.get_genome(params['genome'], self.genome_api)

        self.get_ontology_events(params)
        self.get_translations()
        #
        # # make reports
        summary = self.summarize_gto(params)

        report = self.html_summary(params, summary)
        return report
