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


class CompareAnnotationsUtil:

    workdir = 'tmp/work/'
    staging_dir = "/staging/"
    datadir = "/kb/module/data/"

    translation_locations = {'ec': 'EBI_EC.ModelSEED.json',
                             'keggro': 'KEGG_RXN.ModelSEED.json',
                             'keggko': 'KEGG_KO.ModelSEED.json',
                             'metacyc': 'Metacyc_RXN.ModelSEED.json',
                             'SSO': 'SSO.ModelSEED.json',
                             'go': 'GO.ModelSEED.json'}

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

    def get_genome(self, genome_ref):
        self.genome = self.genome_api.get_genome_v1({"genomes": [{"ref": genome_ref}], 'downgrade': 0})[
            "genomes"][0]['data']

    def get_ontology_events(self, genome_dict):
        if 'ontology_events' in genome_dict:
            for event, ontology in enumerate(genome_dict['ontology_events']):
                self.events[event] = {}
                for term in ontology:
                    self.events[event][term] = ontology[term]
        else:
            logging.info("No ontology events in this genome!")

    def get_translations(self):
        self.translations = {}

        for type in self.translation_locations:
            translations_path = self.datadir + "/" + self.translation_locations[type]
            ontology_translations = json.loads(open(translations_path, "r").read())
            self.translations[type] = {}

            for term in ontology_translations['translation']:
                for entry in ontology_translations['translation'][term]['equiv_terms']:
                    rxn = entry['equiv_term']
                    self.translations[type][term] = rxn

    def search_for_ec(self, line):
        ecList = re.findall(r"\(*[0-9]+\.[0-9\-]+\.[0-9\-]+\.[0-9\-]+", line)
        return(ecList)

    def summarize_gto(self, gto, translations):
        summary = {"genes": {},
                   "terms": {},
                   "rxns": {},
                   "ontology_events": {},
                   "orphan_terms": {}
                   }

        # add ontology events
        for count, oe in enumerate(gto['ontology_events']):
            summary['ontology_events'][count] = oe

        # add gene id to summary
        for feature in gto['features']:
            gene_id = feature['id']
            summary["genes"][gene_id] = {"terms": {},
                                         "rxns": {}
                                         }

            # get ontology term
            if "ontology_terms" in feature:
                for type in feature['ontology_terms']:
                    term_dict = feature['ontology_terms'][type]

                    for term in term_dict:
                        for oe in term_dict[term]:

                            rxn = "none"

                            # get rxn
                            ontology_type = summary['ontology_events'][oe]['id']

                            # fix metacyc terms
                            if ontology_type == 'metacyc':
                                if term.startswith("META:"):
                                    term = term.replace('META:', '')

                            # fix go terms
                            if ontology_type == 'go':
                                if term.startswith("GO:"):
                                    term = term.replace('GO:', '')

                            # fix SSO terms
                            if ontology_type == 'SSO':

                                if term in gto['ontologies_present']['SSO']:
                                    if gto['ontologies_present']['SSO'][term] != 'Unknown':
                                        term = gto['ontologies_present']['SSO'][term]

                            # convert terms to rxns
                            if term in translations[ontology_type]:
                                rxn = translations[ontology_type][term]
                            else:
                                if oe in summary["orphan_terms"]:
                                    summary["orphan_terms"][oe].append(term)
                                    summary["orphan_terms"][oe] = list(
                                        set(summary["orphan_terms"][oe]))
                                else:
                                    summary["orphan_terms"][oe] = [term]

                            # terms
                            if term in summary["genes"][gene_id]['terms']:
                                summary["genes"][gene_id]['terms'][term].append(oe)
                            else:
                                summary["genes"][gene_id]['terms'][term] = [oe]

                            if term in summary['terms']:
                                summary['terms'][term].append(oe)
                                summary['terms'][term] = list(set(summary['terms'][term]))
                            else:
                                summary['terms'][term] = [oe]

                            # rxns
                            if rxn != "none":
                                if rxn in summary["genes"][gene_id]['rxns']:
                                    summary["genes"][gene_id]['rxns'][rxn].append(oe)
                                else:
                                    summary["genes"][gene_id]['rxns'][rxn] = [oe]

                                if rxn in summary['rxns']:
                                    summary['rxns'][rxn].append(oe)
                                    summary['rxns'][rxn] = list(set(summary['rxns'][rxn]))
                                else:
                                    summary['rxns'][rxn] = [oe]

        with open(os.path.join(self.scratch, "summary_dump.json"), 'w') as outfile:
            json.dump(summary, outfile, indent=2)

        return summary

    def html_summary(self, params, summary):

        # convert gto summary for this report
        html_summary_report = {}

        for oe in summary['ontology_events']:
            html_summary_report[oe] = {"gene": [], "term": [], "rxn": []}

        for gene in summary["genes"]:
            for term in summary["genes"][gene]['terms']:
                for oe in summary["genes"][gene]['terms'][term]:
                    html_summary_report[oe]['gene'].append(gene)
                    html_summary_report[oe]['term'].append(term)

                    html_summary_report[oe]['gene'] = list(set(html_summary_report[oe]['gene']))
                    html_summary_report[oe]['term'] = list(set(html_summary_report[oe]['term']))

            for rxn in summary["genes"][gene]['rxns']:
                for oe in summary["genes"][gene]['rxns'][rxn]:
                    html_summary_report[oe]['rxn'].append(rxn)
                    html_summary_report[oe]['gene'].append(gene)

                    html_summary_report[oe]['rxn'] = list(set(html_summary_report[oe]['rxn']))
                    html_summary_report[oe]['gene'] = list(set(html_summary_report[oe]['gene']))

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

    def run(self, ctx, params):

        logging.info(params['annotations_to_compare'])

        # collect some metadata
        self.get_genome(params['genome'])
        self.get_ontology_events(self.genome)
        self.get_translations()

        # make reports
        summary = self.summarize_gto(self.genome, self.translations)

        report = self.html_summary(params, summary)
        return report

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
                'description', summary["ontology_events"][o]['method']) + '_' + str(o)

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
