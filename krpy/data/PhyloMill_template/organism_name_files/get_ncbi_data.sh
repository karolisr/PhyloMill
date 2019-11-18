#!/usr/bin/env bash

curl -s -o taxdump.tar.gz ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
mkdir ncbi
tar xzf taxdump.tar.gz -C ncbi
rm -rf taxdump.tar.gz
mv ncbi/names.dmp ncbi_tax_names
rm -rf ncbi
