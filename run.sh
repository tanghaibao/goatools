#!/bin/bash

#python obo_parser.py gene_ontology.1_2.obo
python genemerge.py --alpha 0.05 --hierarchy data/association data/description data/accns data/alpha_gtp95.list stdout 
