import pandas as pd
import sys
import json

annotations_file_path = "/Users/kimbrel1/Github/MergeMetabolicAnnotations/test/test_data/evidence_test.txt"

params = 0


##
annotations = annotations.drop_duplicates()

ontology = {
    'event_id': 'description',
    'description': 'description',
    'ontology_id': 'ontology',
    'ontology_terms': {}
}

# add imported terms
for index, row in annotations.iterrows():
    if pd.notnull(row['term']):
        if row['gene'] in ontology['ontology_terms']:
            if params == 0:
                ontology['ontology_terms'][row['gene']].append(
                    {'term': row['term']}
                )
            if params == 1:
                evidence_terms = row.to_dict()
                evidence_terms.pop('gene', None)
                evidence_terms.pop('term', None)

                ontology['ontology_terms'][row['gene']].append(
                        {'term': row['term'],
                        'evidence': {'scores': evidence_terms}}
                )
        else:
            if params == 0:
                ontology['ontology_terms'][row['gene']] = [
                        {'term': row['term']}
                    ]
            if params == 1:

                evidence_terms = row.to_dict()
                evidence_terms.pop('gene', None)
                evidence_terms.pop('term', None)

                ontology['ontology_terms'][row['gene']] = [
                    {'term': row['term'],
                        'evidence': {'scores': evidence_terms}}
                ]

                

                
                
                
                # print(json.dumps(pivoted, indent=2))
                


    






print(json.dumps(ontology, indent=2))