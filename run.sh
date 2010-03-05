#!/bin/bash

echo "select the test from below"
select TEST in `seq 1 3`;
do
case $TEST in

1)
python scripts/find_enrichment.py --alpha=0.05 --indent data/study data/population data/association
;;

2)
python scripts/plot_go_term.py --term=GO:0008135
;;

3)
python scripts/plot_go_network.py data/go_counts
;;

esac
done
