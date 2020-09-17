import sys
import json
import pandas as pd

translations = {"translation": {},
                "comment": "SSO->MS TranslationTable",
                "ontology1": "sso"}
tsv = pd.read_csv(sys.argv[1], sep="\t", names=['term', 'translation'])

for index, row in tsv.iterrows():
    if row['term'] not in translations['translation']:
        translations['translation'][row['term']] = {"equiv_terms": []}

    translations['translation'][row['term']]["equiv_terms"].append(
        {"equiv_term": row['translation']})


print(json.dumps(translations, indent=4, sort_keys=True))
