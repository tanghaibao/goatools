#!/bin/bash

# configuration for Gene Ontology files
GO_OBO_FILE=go-basic.obo
GOSLIM_OBO_FILE=goslim_generic.obo

GO_OBO_DOWNLOAD=http://purl.obolibrary.org/obo/go/go-basic.obo
GOSLIM_OBO_DOWNLOAD=http://www.geneontology.org/ontology/subsets/goslim_generic.obo

# if the gene ontology files don't exist, download them
if [ ! -f $GO_OBO_FILE ]
then
    echo "downloading GO file: $GO_OBO_FILE"
    wget -O $GO_OBO_FILE $GO_OBO_DOWNLOAD
fi

if [ ! -f $GOSLIM_OBO_FILE ]
then
    echo "downlaoding GOslim file: $GOSLIM_OBO_FILE"
    wget -O $GOSLIM_OBO_FILE $GOSLIM_OBO_DOWNLOAD
fi

echo "select the test from below"
select TEST in 'Test enrichment' 'Test plotting go terms' 'Test mapslim' 'Test mapslim on assocation'


do
case $REPLY in

1)
python scripts/find_enrichment.py --alpha=0.05 --indent data/study data/population data/association
;;

2)
python scripts/plot_go_term.py --term=GO:0008135
;;

3)
python 'tests/test_mapslim.py'
;;

4)
python scripts/map_to_slim.py --association_file=data/association --slim_out=direct $GO_OBO_FILE $GOSLIM_OBO_FILE
;;


esac
done
