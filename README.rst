Tools for Gene Ontology
========================

:Author: Haibao Tang (tanghaibao), Brent Pedersen (brentp), Aurélien Naldi
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
- Python version >= 2.6, try install this package first. Within this folder::

    easy_install .

- ``.obo`` file for the most current `gene ontology <http://www.geneontology.org/>`_::

    wget http://geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo 

- `fisher <http://pypi.python.org/pypi/fisher/>`_ module for calculating Fisher's exact test::
    
    easy_install fisher

And put it in the current folder where you want to do your analysis.

If you need to plot the ontology lineage, you need the following to tools to be installed.

- `Graphviz <http://www.graphviz.org/>`_, for graph visualization.

- `pygraphviz <http://networkx.lanl.gov/pygraphviz/>`_, Python binding for communicating with Graphviz::

    easy_install pygraphviz 


Cookbook
---------
``run.sh`` contains example cases, which calls the utility scripts in the ``scripts`` folder.

Find GO enrichment of genes under study
::::::::::::::::::::::::::::::::::::::::::
see ``find_enrichment.py`` for usage. It takes as arguments files containing:

* gene names in a study
 
* gene names in population (or other study if --compare is specified)

* an association file that maps a gene name to a GO category.

please look at ``data/`` folder to see examples on how to make these files. when ready, the command looks like::

    python scripts/find_enrichment.py --alpha=0.05 --indent data/study data/population data/association

and can filter on the significance of erichment or purification.
it can report various multiple testing corrected p-values as well as
the false discovery rate.

Read and plot GO lineage
::::::::::::::::::::::::::::::::::::
see ``plot_go_term.py`` for usage. 

``plot_go_term.py`` can plot the lineage of a certain GO term, by::

   python scripts/plot_go_term.py --term=GO:0008135

this will plot the following image.

.. image:: http://lh6.ggpht.com/_srvRoIok9Xs/S9HhleQrk5I/AAAAAAAAA5U/dzVIvjlYCQU/s800/GO_0008135.png 
    :alt: term lineage


