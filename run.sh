#!/bin/bash

#python goatools/obo_parser.py gene_ontology.1_2.obo -t GO:0008135
python goatools/go_enrichment.py --alpha 0.05 --indent data/study data/population data/association
