import json
import numpy as np
from tqdm import tqdm
import openai
import requests
import json
import tempfile
import urllib.request
import shutil
import torch
from tqdm import tqdm
from xml.etree import ElementTree
from transformers import BertModel, BertTokenizer
from google.cloud import storage
import pickle
from openai import OpenAI
import time

from itertools import islice
# remember to set your OpenAI api key
client = OpenAI(
  api_key='api_key',
)

def get_embedding(text, model="text-embedding-ada-002"):
    text = text.replace("\n", " ")
    try:
        output = client.embeddings.create(input = [text], model=model).data[0].embedding
    except Exception as e:
        output = np.zeros(1536)
        print(e)
        print("##################")
        print(text)
        print("##################")
    return output


def read_json_url(url: str) -> dict:
    """
    Read a JSON-formatted file from a URL.

    Example:
        >>> data = read_json_url(
        ...     "https://storage.googleapis.com/my-bucket/my-data.json"
        ... )
        >>> print(data)  # The JSON content as a dictionary

    Args:
        url: URL of the JSON file.

    Returns:
        A dictionary representation of the JSON file.
    """
    # Ensure the URL uses a valid scheme
    assert url.startswith(('http://', 'https://', 'ftp://')), 'URL should start with http://, https://, or ftp://'
    
    # Open the URL and create a temporary file to hold its contents
    with urllib.request.urlopen(url) as response:
        with tempfile.NamedTemporaryFile(mode='wb+', delete=False) as tmp_file:  # Change 'w+' to 'wb+'
            shutil.copyfileobj(response, tmp_file)  # Copy the content of the response to the temporary file
            tmp_file.seek(0)  # Go to the start of the file before reading
            return json.load(tmp_file)  # Load the JSON content as a dictionary

Summary_path = "data/ncbi_summary_documents.json"
Gene_vocab_path = "data/biomart_gene2ensembl.json"
#NCBI_cleaned_summary_of_genes = read_json_url(file_path)
with open(Gene_vocab_path, 'r') as json_file:
            data = json.load(json_file)
            Gene_IDs = data.keys()

gpt_gene_name_to_embedding = {}


GPT_DIM = 1536 # fix GPT embeddings
for key in tqdm(Gene_IDs):
    if key not in gpt_gene_name_to_embedding:
        #print('key',key)
        if Gene_IDs[key] == '': 
            # if the dictionary does not have information about a gene
            gpt_gene_name_to_embedding[key] = np.zeros(GPT_DIM) # it's hard coded to be 0
        else:
            text = (
                f"what's the functional annotation, molecular mechanism, "
                f"protein characteristic, protein-protein interaction, and RNA studies, "
                f"epigenetic studies, GWAS study results about gene {Gene_IDs[key]}"
            )  # this is the question
            if len(text.split()) > 4000:
                text = " ".join(text.split()[:4000])
            #print(text)
            if not text:
                # If the dictionary does not have information about the gene
                embedding = np.zeros(GPT_DIM)
            else:
                embedding = get_embedding(text)
                gpt_gene_name_to_embedding[key] = embedding
                #print(embedding)
        
    # save NCBI_cleaned_summary_of_genes
    with open(f"data/NCBI_summary_OpenAI_Embedding.pickle", "wb") as fp:
        pickle.dump(gpt_gene_name_to_embedding_clean_text, fp)
        # sleep for 1 second
    time.sleep(1)

