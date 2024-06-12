#!/usr/bin/env python3

###############################################################################
# script name: get_pubmed_attributes.py
# developed by: Savvas Paragkamian
# framework: Crete IDEA
###############################################################################
# GOAL:
# Aim of this script is to retrieve information from pubmed about mesh terms
# and year of publication based on PMIDs and output them as tsv
###############################################################################
# usage:./scripts/get_pubmed_attributes.py file/with/PMIDs path/to/output_file.tsv
###############################################################################
import xml.etree.ElementTree as ET
import requests
import csv
import os,sys, time

###################### user input ##################
pmids_files = sys.argv[1]
output_file = sys.argv[2]

pmids = []

with open(pmids_files) as file:
    for line in file:
        
        l=line.split('\t')
        pmid = l[0].replace("PMID:","")
        pmids.append(pmid)

# make unique
pmids = list(set(pmids))
print("there are " + str(len(pmids)) + " pmids")

#################### request script ################
def attr_request(accession):
    """
    Retrieves titles and abstracts for a list of PubMed identifiers.
    """
    
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    
    params = {
        'db': 'pubmed',
        'id': accession,
        'retmode': 'xml'}
    xml_data = {}
    mesh_terms = []
    
    start_time = time.time()
    connection_timeout = 60 # seconds

    while True:
        try:
            response = requests.get(base_url, params=params)
            break
        except requests.ConnectionError:
            if time.time() > start_time + connection_timeout:
                raise Exception('Unable to get updates after {} seconds of ConnectionErrors'.format(connection_timeout))
            else:
                time.sleep(1)

    code = response.status_code

    if code == 204:
        print(str(accession) + "\t" + str(code))

    if code == 400:
        print(str(accession) + "\t" + str(code))

    if code == 200:
        print(str(accession) + "\t" + str(code))
        article = ET.fromstring(response.content)
    
        xml_data['pmid'] = article.findtext('.//PMID')
        xml_data['title'] = article.findtext('.//ArticleTitle')
        xml_data['abstract'] = article.findtext('.//AbstractText')
        xml_data['year'] = article.findtext('.//Year')
        for mesh_heading in article.findall(".//MeshHeading"):
            descriptor_name = mesh_heading.find("DescriptorName")
            if descriptor_name is not None:
                mesh_terms.append(descriptor_name.text)
            else:
                mesh_terms_str = "NA"
                # Join MeSH terms into a single string
        xml_data['mesh_terms'] = "; ".join(mesh_terms)
            
            
        # Write the dictionary to the data list
        return(xml_data)

    time.sleep(uniform(1,2))

        
################################## main code ################################
articles = []

#Without an API key, any site (IP address) posting more than
#3 requests per second to the E-utilities will receive an error message.
#By including an API key, a site can post up to 10 requests per second by default.

fieldnames = ['pmid','title', 'abstract', 'year', 'mesh_terms']

with open(output_file, 'w', newline='') as tsvfile:
    writer = csv.DictWriter(tsvfile,fieldnames=fieldnames,delimiter='\t')
    
    # Write the header
    writer.writeheader() 
    # Write the data rows per dictionary/file

    for ref in pmids:
        xml_data = attr_request(ref)

        writer.writerow(xml_data)

    
print(f"Conversion complete. The TSV file '{output_file}' has been created.")

