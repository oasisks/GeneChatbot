from bs4 import BeautifulSoup
import requests
import pprint
import re


def _gene_parser(gene_html):
    response = requests.get(gene_html)
    gene_soup = BeautifulSoup(response.text, features="html.parser")

    summary = gene_soup.find("dl", {"id": "summaryDl"})

    dts = summary.findAll("dt")[:-2]
    dds = summary.findAll("dd")[:-2]

    # there needs to be equal amounts of title and description
    if len(dts) != len(dds):
        return None

    result = {}

    for i in range(len(dts)):
        dt = dts[i]
        dd = dds[i]

        dt_re = None
        dd_re = None
        for dt_content in dt.contents:
            if dt_content != "\n":
                res = re.findall(r'\w+', dt_content)
                dt_re = " ".join(res)
                break

        for dd_content in dd.contents:
            if dd_content != "\n":
                if not isinstance(dd_content, str):
                    dd_content = dd_content.contents[0]
                res = re.findall(r'\w+(?:-\w+)*', dd_content)
                dd_re = " ".join(res)
                break

        if dt_re is not None and dd_re is not None:
            result[dt_re] = dd_re

    result["url"] = gene_html
    return result


def get_gene_data(gene_name: str):
    """
    Given the gene_name, it will return the summary of the gene
    :param gene_name: the name of the gene symbol
    :return:
    """
    url = f"https://www.ncbi.nlm.nih.gov/gene/?term={' '.join(gene_name.split())}"
    gene_id_url = "https://www.ncbi.nlm.nih.gov/gene/"
    response = requests.get(url)
    soup = BeautifulSoup(response.text, features="html.parser")
    gene_name_ids = soup.findAll("td", {"class": "gene-name-id"})

    genes = {}

    for gene_name_id in gene_name_ids:
        gene_id = gene_name_id.find("span", {"class": "gene-id"}).contents[0].strip("ID: ")
        result = _gene_parser(gene_id_url + gene_id)
        genes[gene_id] = result

    return genes


def get_disease_paper(gene_name: str):
    """
    Given the gene_name, it will try to query the website for information and return an appropriate response
    The returned result will most likely be the first query on the search query.
    :param gene_name:
    :return:
    """
    pass


if __name__ == '__main__':
    result = get_gene_data("PRKN")
    pprint.pprint(result)
