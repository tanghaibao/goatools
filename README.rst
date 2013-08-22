Tools for Gene Ontology
========================

:Author: Haibao Tang (tanghaibao), Brent Pedersen (brentp), Aurelien Naldi (aurelien-naldi)
:Email: tanghaibao@gmail.com
:License: BSD

.. contents ::

Description
------------
This package contains a Python library to

- process over- and under-representation of certain GO terms, based on Fisher's exact test. Also implemented several multiple correction routines (including Bonferroni, Sidak, and false discovery rate).
- process the obo-formatted file from `Gene Ontology website <http://geneontology.org>`_. The data structure is a directed acyclic graph (DAG) that allows easy traversal from leaf to root.
- map GO terms (or protein products with multiple associations to GO terms) to GOslim terms (analog to the map2slim.pl script supplied by geneontology.org)


Installation
-------------
- Python version >= 2.6, try install this package first. Within this folder::

    easy_install .

- ``.obo`` file for the most current `gene ontology <http://www.geneontology.org/>`_::

    wget http://geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo

- ``.obo`` file for the most current `GO Slim <http://www.geneontology.org/GO.slims.shtml>`_ terms (.e.g generic GOslim) ::

    wget http://www.geneontology.org/GO_slims/goslim_generic.obo

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

please look at ``tests/data/`` folder to see examples on how to make these files. when ready, the command looks like::

    python scripts/find_enrichment.py --pval=0.05 --indent data/study data/population data/association

and can filter on the significance of (e)nrichment or (p)urification.
it can report various multiple testing corrected p-values as well as
the false discovery rate.

The "e" in the "Enrichment" column means "enriched" - the concentration of GO term in
the study group is significantly *higher* than those in the population.  The "p" stands
for "purified" - significantly *lower* concentration of the GO term in the study group
than in the population.


Read and plot GO lineage
::::::::::::::::::::::::::::::::::::
see ``plot_go_term.py`` for usage.

``plot_go_term.py`` can plot the lineage of a certain GO term, by::

   python scripts/plot_go_term.py --term=GO:0008135

this will plot the following image.

.. image:: http://lh6.ggpht.com/_srvRoIok9Xs/S9HhleQrk5I/AAAAAAAAA5U/dzVIvjlYCQU/s800/GO_0008135.png
    :alt: GO term lineage

Sometimes people like to stylize the graph themselves, use option ``--gml`` to
generate a GML output which can then be used in an external graph editing
software like `Cytoscape <http://www.cytoscape.org/>`_. The following image is
produced by importing the GML file into Cytoscape using yFile orthogonal
layout and solid VizMapping. Note that the `GML reader plugin
<https://code.google.com/p/graphmlreader/>`_ may need to be
downloaded and installed in the ``plugins`` folder of Cytoscape::

    python scripts/plot_go_term.py --term=GO:0008135 --gml

.. image:: http://tinyurl.com/by2m57n
    :alt: GO term lineage (Cytoscape)


Map GO terms to GOslim terms
::::::::::::::::::::::::::::::::::::
see ``map_to_slim.py`` fro usage. As arguments it takes the gene ontology files:

* the current gene ontology file ``gene_ontology.1_2.obo``

* the GOslim file to be used (e.g. ``goslim_generic.obo`` or any other GOslim file)

The script either maps one GO term to it's GOslim terms, or protein products with multiple associations to all it's GOslim terms.

To determine the GOslim terms for a single GO term, you can use the following
command::

    python scripts/map_to_slim.py --term=GO:0008135 gene_ontology.1_2.obo goslim_generic.obo

To determine the GOslim terms for protein products with multiple associations::

    python scripts/map_to_slim.py --association_file=data/association gene_ontology.1_2.obo goslim_generic.obo

Where the ``association`` file has the same format as used for ``find_enrichment.py``.

The implemented algorithm is described in more detail at the go-perl documenation of `map2slim
<http://search.cpan.org/~cmungall/go-perl/scripts/map2slim>`_.