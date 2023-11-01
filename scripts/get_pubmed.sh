#!/bin/bash

cd ~/pubmed2023

wget -m -nH --cut-dirs=2 'ftp://ftp.ncbi.nlm.nih.gov/pubmed/baseline/'
wget -m -nH --cut-dirs=2 'ftp://ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/'

for i in `ls -1 *.xml.gz`; do
    echo $i
    suffix='.tsv.gz'
    FILE=${i/.xml.gz/$suffix} 
    tsv=${FILE/.tsv.gz/'.tsv'}

    if [ ! -f $FILE ];then
        gzip -cd $i | /data/databases/scripts/gathering_data/pubmed/pubmed.pl > output_temp.txt
        gzip -9 output_temp.txt
        mv output_temp.txt.gz $FILE
        echo $FILE " is ready"
    fi
done
