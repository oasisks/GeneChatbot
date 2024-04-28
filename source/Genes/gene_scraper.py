import Bio.Entrez.Parser
from bs4 import BeautifulSoup
from dotenv import load_dotenv
from pymongo.mongo_client import MongoClient
import os
from Bio import Entrez
import json
from tqdm import tqdm
import pprint

load_dotenv()
NCBI_KEY = os.getenv("NCBI_API_KEY")
MONGO_API_KEY = os.getenv("MONGO_API_KEY")
Entrez.email = "something@something.com"
Entrez.api_key = NCBI_KEY
uri = os.getenv("MONGO_URI")


def get_ncbi_gene_id(ensemble_id: str, name: str):
    """
    Returns the gene_id in the NCBI website

    If there are multiple IDs, then it will return the first ID
    :param ensemble_id: the universal id of the gene
    :param name: the name associated with the ensemble_id
    :return: a list of ids
    """
    handle = Entrez.esearch(db="gene", term=ensemble_id)
    record = Entrez.read(handle)
    handle.close()
    ids = record["IdList"]
    handle = Entrez.esearch(db="gene", term=name)
    handle.close()
    ids = record["IdList"]
    return ids


def get_ncbi_gene_data(gene_id: Bio.Entrez.Parser.StringElement):
    """
    Grabs all the information on the NCBI website given the NCBI Gene Id
    :param gene_id: the id of the gene
    :return:
    """
    handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
    record = handle.read()
    handle.close()
    return record


def ncbi_gene_parser(markdown: bytes):
    """
    This returns a dictionary of entries related to the gene
    :param markdown: a xml markdown of the gene data
    :return: Returns a list of information related to the gene
    """
    symbol_report_url = "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:"
    soup = BeautifulSoup(markdown, features="xml")

    # unique tags
    lineages = soup.find("OrgName_lineage")
    genus = soup.find("BinomialOrgName_genus")
    species = soup.find("BinomialOrgName_species")
    gene_type = soup.find("Entrezgene_type")

    # under references
    gene_ref = soup.find("Gene-ref")
    official_symbol = gene_ref.find("Gene-ref_locus")
    official_full_name = gene_ref.find("Gene-ref_desc")

    # also known as
    ref_syns = soup.findAll("Gene-ref_syn_E")
    similar_names = ref_syns if ref_syns is None else [elt.text for elt in ref_syns]

    # summary
    summary = soup.find("Entrezgene_summary")

    # GeneRIF
    gene_rifs = soup.findAll("Gene-commentary")
    title_to_id = []
    pubmed_ids = set()
    if gene_rifs is not None:
        for gene_commentary in gene_rifs:
            title = gene_commentary.find("Gene-commentary_text")
            pubmed_id = gene_commentary.find("PubMedId")
            if title is not None and pubmed_id is not None and pubmed_id.text not in pubmed_ids:
                title_to_id.append({
                    "title": title.text,
                    "pubMedId": int(pubmed_id.text)
                })
                pubmed_ids.add(pubmed_id.text)

    return {
        "organism": {
            "genus": genus if genus is None else genus.text,
            "species": species if species is None else species.text
        },
        "lineages": lineages if lineages is None else lineages.text.split(";"),
        "reference": gene_ref if gene_ref is None else {
            "symbol": official_symbol if official_symbol is None else official_symbol.text,
            "official_name": official_full_name if official_full_name is None else official_full_name.text
        },
        "similar_names": None if not similar_names else similar_names,
        "summary": summary if summary is None else summary.text,
        "gene_type": gene_type if gene_type is None else gene_type["value"],
        "gene_rifs": gene_rifs if gene_rifs is None else title_to_id
    }


def main():
    file = open("../Data/biomart_ensembl2gene.json")
    document = json.load(file)
    file.close()
    mongo_client = MongoClient(uri)

    try:
        mongo_client.admin.command('ping')
        print("Pinged your deployment. You successfully connected to MongoDB!")
    except Exception as e:
        print(e)
    gene_db = mongo_client["gene_db"]
    raw_data = gene_db["raw"]
    ids = set()
    get_ncbi_gene_id("dasd", "5.8S-rRNA")
    for ensemble_id, name in document.items():
        # count = raw_data.find({"_id": ensemble_id}).collection.count_documents({})
        # if count > 0:
        #     continue
        # gene_id = get_ncbi_gene_id(ensemble_id, name)
        # if there are no gene ids or the id is already queried, we continue
        pass
        # if not len(gene_id) > 0 or gene_id[0] in ids:
        #     continue
        #
        # markdown = get_ncbi_gene_data(gene_id[0])
        # data = ncbi_gene_parser(markdown)
        # data["_id"] = ensemble_id
        # response = raw_data.insert_one(data)
        # ids.add(gene_id[0])


if __name__ == '__main__':
    main()

