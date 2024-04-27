import json
import time
from Bio import Entrez
from xml.etree import ElementTree as ET
from tqdm import tqdm
import os
from requests import HTTPError

# Entrez email and API key
Entrez.email = "something@something.com"
Entrez.api_key = open("api_key").read().strip()

def entrez_search(query):
    handle = Entrez.esearch(db="gene", term=query)
    record = Entrez.read(handle)
    handle.close()
    return record

def extract_gene_ids(esearch_result):
    gene_ids = esearch_result["IdList"]
    assert len(gene_ids) == 1, "more than one gene id found"
    return gene_ids[0]


def fetch_gene_information(gene_id):
    try:
        handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
        record = handle.read()
        handle.close()
    except HTTPError as http_err:
        if http_err.response.status_code == 429 or http_err.response.status_code == 400:
            print("Too many requests. Sleeping for 10 seconds.")
            time.sleep(10)
            return fetch_gene_information(gene_id)
        else:
            raise


    summary, commentaries_text = extract_gene_summary_and_commentaries(record)
    return summary, commentaries_text

def extract_gene_summary_and_commentaries(efetch_result):
    root = ET.fromstring(efetch_result)
    gene_summary = root.find(".//Entrezgene_summary")
    summary = gene_summary.text if gene_summary is not None else None
    generif_commentary_texts = [elem.find("Gene-commentary_text").text
                                for elem in root.findall(".//Gene-commentary")
                                if elem.find("Gene-commentary_text")!=None and
                                elem.find("Gene-commentary_type").text == '18']
    generif_commentary_texts = [text for text in generif_commentary_texts if len(text) > 50]
    return summary, generif_commentary_texts


#Load the ensembl dictionary which includes all the genes with their ENSG equivalent
with open("data/vocab/biomart_gene2ensembl.json") as json_file:
    gene_vocab = json.load(json_file)
    gene_vocab = list(set(gene_vocab.values()))


base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

'''
{
    'ensemble_id': {
        'gene_name': 'A1BG',
        'summary: 'dahfdsah',
        'commentaries': ['dahfdsah', 'dahfdsah']
    }

}

'''

if os.path.exists("data/ncbi_documents.json"):
    with open("data/ncbi_documents.json") as json_file:
        current_ncbi_documents = json.load(json_file)
else:
    current_ncbi_documents = {}

# Iterate through the gene vocabulary to fetch their summaries and commentaries
gene_data_ncbi = current_ncbi_documents
i = 0
gene_vocab = sorted(gene_vocab, reverse=False)
for gene_name in tqdm(gene_vocab, total=len(gene_vocab)):
    # if i>50000:
    #     break
    i+=1


    ensemble_id = gene_name
    if ensemble_id in current_ncbi_documents:
        continue


    try:
        esearch_result = entrez_search(ensemble_id)
        gene_id = extract_gene_ids(esearch_result)

        commentaries = []
        summary = None
        if gene_id:
            summary, commentaries = fetch_gene_information(gene_id)

        gene_data_ncbi[ensemble_id] = {
            'ensemble_id': gene_name,
            'summary': summary,
            'commentaries': commentaries
        }

        with open("data/ncbi_documents.json", "w") as json_file:
            json.dump(gene_data_ncbi, json_file, indent=4)

    except Exception as e:
        print("An error occurred:", str(e))

    time.sleep(1)


