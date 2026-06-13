import json


def get_event_lists(ontology):
    # print(type(ontology['feature_types']))

    #gene_features = [k for k, v in ontology['feature_types'].items() if v == "gene"]

    events = {}
    for event in ontology["ontology_events"]:
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

        # print(events)

        for gene in event["ontology_terms"]:
            print(gene)
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


ontology_selected = json.loads(
        open("/Users/kimbrel1/Library/CloudStorage/Dropbox/Downloads/LLNL/Cb_RPG.JSON/9.json", "r").read())


get_event_lists(ontology_selected)