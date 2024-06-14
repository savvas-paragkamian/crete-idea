#!/usr/bin/env python3

###############################################################################
# script name: pubmed_mesh_xml_to_tsv.py
# developed by: Savvas Paragkamian
# framework: Crete IDEA
###############################################################################
# GOAL:
# Aim of this script is to transform all xml structure of mesh terms 
# of pubmed to xml. There are about 30000 terms
###############################################################################
# usage:./scripts/pubmed_mesh_xml_to_tsv.py 
###############################################################################

import xml.etree.ElementTree as ET
import csv
import os,sys, time

root = ET.parse("data/desc2024.xml")
# the mesh terms xml was downloaded from https://www.nlm.nih.gov/databases/download/mesh.html
# the mesh xml has DescriptorRecord as the begining of each record

fields = ['DescriptorUI', 'DescriptorName', 'TreeNumber']

with open('results/pubmed_mesh_terms.tsv', 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t')
    
    # Write the header
    writer.writerow(fields)
    # parsing the xml file

    for record in root.findall('DescriptorRecord'):
        # Extract descriptors
        descriptor_ui = record.find('DescriptorUI').text
        descriptor_name = record.find('DescriptorName/String').text
    
        # Extract tree numbers
        tree_numbers = [tree_number.text for tree_number in record.findall('TreeNumberList/TreeNumber')]
        tree_numbers_str = "; ".join(tree_numbers)
        
        # Write to file
        writer.writerow([descriptor_ui,descriptor_name,tree_numbers_str])


