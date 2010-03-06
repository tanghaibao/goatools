Tools for Gene Ontology
========================

:Author: Haibao Tang (tanghaibao), Brent Pedersen (brentp)
:Email: tanghaibao@gmail.com
:License: BSD

.. contents ::

Description
------------
This package contains a Python library to

- process over- and under-representation of certain GO terms, based on Fisher's exact test. Also implemented several multiple correction routines (including Bonferroni, Sidak, and false discovery rate).
- process the obo-formatted file from `Gene Ontology website <http://geneontology.org>`_. The data structure is a directed acyclic graph (DAG) that allows easy traversal from leaf to root.


Installation
-------------
- ``.obo`` file for the most current `gene ontology <http://www.geneontology.org/>`_::

    wget http://geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo 

And put it in the current folder where you want to do your analysis.

If you need to plot the ontology lineage, you need the following to tools to be installed.

- `Graphviz <http://www.graphviz.org/>`_, for graph visualization.

- `pygraphviz <http://networkx.lanl.gov/pygraphviz/>`_, Python binding for communicating with Graphviz::

    easy_install pygraphviz 


Cookbook
---------
``run.sh`` contains a few example usages, which calls the utility scripts in the ``scripts`` folder.

Find GO enrichment of genes under study
::::::::::::::::::::::::::::::::::::::::::
see ``find_enrichment.py`` for usage. 


Read and plot GO lineage
::::::::::::::::::::::::::::::::::::
see ``plot_go_term.py`` and ``plot_go_network.py`` for usage. 

``plot_ge_term.py`` can plot the lineage of a certain GO term, by::

   python scripts/plot_go_term.py --term=GO:0008135

this will plot the following image.

.. image:: http://img502.imageshack.us/img502/7016/go0008135.png 
    :alt: term lineage


