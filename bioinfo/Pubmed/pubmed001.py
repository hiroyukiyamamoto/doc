# coding: utf-8
import urllib.request
from xml.etree.ElementTree import *

keyword = 'metabolomics AND "machine learning"'
baseURL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term="


#https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=metabolomics%20AND%20%22machine%20learning%22

def get_id(url):#論文IDを取得する
    result = urllib.request.urlopen(url)
    return result

def main():
    url = baseURL + keyword
    result = get_id(url)
    element = fromstring(result.read())
    filename = "idlist_"+keyword+".txt"
    f = open(filename, "w")
    for e in element.findall(".//Id"):
        f.write(e.text)
        f.write("\n")
    f.close()

if __name__ == "__main__":
    main()