Tools for Gene Ontology
========================

.. image:: https://pypip.in/v/goatools/badge.png
    :target: https://crate.io/packages/goatools/
    :alt: Latest PyPI version

.. image:: https://pypip.in/d/goatools/badge.png
    :target: https://crate.io/packages/goatools/
    :alt: Number of PyPI downloads

.. image:: https://pypip.in/license/goatools/badge.png
    :target: https://crate.io/packages/goatools/
    :alt: License

.. image:: https://scrutinizer-ci.com/g/tanghaibao/goatools/badges/quality-score.png
    :target: https://scrutinizer-ci.com/g/tanghaibao/goatools
    :alt: Scrutinizer

.. image:: https://scrutinizer-ci.com/g/tanghaibao/jcvi/badges/build.png
    :target: https://scrutinizer-ci.com/g/tanghaibao/goatools
    :alt: Build

:Author: Haibao Tang (`tanghaibao <http://github.com/tanghaibao>`_),
         Brent Pedersen (`brentp <http://github.com/brentp>`_),
         Aurelien Naldi (`aurelien-naldi <http://github.com/aurelien-naldi>`_),
         Patrick Flick (`r4d2 <http://github.com/r4d2>`_),
         Jeff Yunes (`yunesj <http://github.com/yunesj>`_),
         Kenta Sato (`bicycle1885 <http://github.com/bicycle1885>`_),
         Chris Mungall (`cmungall <https://github.com/cmungall>`_),
         Greg Stupp (`stuppie <https://github.com/stuppie>`_),
         Debra Klopfenstein(`dvklopfenstein <https://github.com/dvklopfenstein>`_)
:Email: tanghaibao@gmail.com
:License: BSD

.. contents ::

Description
------------
This package contains a Python library to

- process over- and under-representation of certain GO terms, based on Fisher's
  exact test. Also implemented several multiple correction routines (including
  Bonferroni, Sidak, and false discovery rate).
- process the obo-formatted file from `Gene Ontology website <http://geneontology.org>`_.
  The data structure is a directed acyclic graph (DAG) that allows easy traversal
  from leaf to root.
- map GO terms (or protein products with multiple associations to GO terms) to
  GOslim terms (analog to the map2slim.pl script supplied by geneontology.org)


Installation
-------------
- Make sure your Python version >= 2.6, install it via PyPI::

    easy_install goatools

- ``.obo`` file for the most current `GO <http://geneontology.org/page/download-ontology>`_::

    wget http://purl.obolibrary.org/obo/go/go-basic.obo

- ``.obo`` file for the most current `GO Slim <http://geneontology.org/page/go-slim-and-subset-guide>`_
  terms (.e.g generic GOslim) ::

    wget http://www.geneontology.org/ontology/subsets/goslim_generic.obo


Dependencies
-------------
- `fisher <http://pypi.python.org/pypi/fisher/>`_ module for calculating
  Fisher's exact test::

    easy_install fisher

- To plot the ontology lineage, install one of these two options:

  1. Graphviz

     - `Graphviz <http://www.graphviz.org/>`_, for graph visualization.
     - `pygraphviz <http://networkx.lanl.gov/pygraphviz/>`_, Python binding for
       communicating with Graphviz::

         easy_install pygraphviz

  2. `pydot <https://code.google.com/p/pydot/>`_, a Python interface to Graphviz's Dot language.

     * `pyparsing <http://pyparsing.wikispaces.com/>`_ is a prerequisite for pydot
     * Images can be viewed using either:

       * `ImageMagick <http://www.imagemagick.org/>`_'s *display*
       * `Graphviz <http://www.graphviz.org/>`_


Cookbook
---------
``run.sh`` contains example cases, which calls the utility scripts in the
``scripts`` folder.

Find GO enrichment of genes under study
::::::::::::::::::::::::::::::::::::::::::
See ``find_enrichment.py`` for usage. It takes as arguments files containing:

* gene names in a study
* gene names in population (or other study if --compare is specified)
* an association file that maps a gene name to a GO category.

Please look at ``tests/data/`` folder to see examples on how to make these
files. when ready, the command looks like::

    python scripts/find_enrichment.py --pval=0.05 --indent data/study data/population data/association

and can filter on the significance of (e)nrichment or (p)urification.
it can report various multiple testing corrected p-values as well as
the false discovery rate.

The "e" in the "Enrichment" column means "enriched" - the concentration of GO
term in the study group is significantly *higher* than those in the population.
The "p" stands for "purified" - significantly *lower* concentration of the GO
term in the study group than in the population.


Read and plot GO lineage
::::::::::::::::::::::::::::::::::::
See ``plot_go_term.py`` for usage.  ``plot_go_term.py`` can plot the lineage of
a certain GO term, by::

   python scripts/plot_go_term.py --term=GO:0008135

This command will plot the following image.

.. image:: https://dl.dropboxusercontent.com/u/15937715/Data/github/goatools/gograph.png
    :alt: GO term lineage

Sometimes people like to stylize the graph themselves, use option ``--gml`` to
generate a GML output which can then be used in an external graph editing
software like `Cytoscape <http://www.cytoscape.org/>`_. The following image is
produced by importing the GML file into Cytoscape using yFile orthogonal
layout and solid VizMapping. Note that the `GML reader plugin
<https://code.google.com/p/graphmlreader/>`_ may need to be
downloaded and installed in the ``plugins`` folder of Cytoscape::

    python scripts/plot_go_term.py --term=GO:0008135 --gml

.. image:: https://dl.dropboxusercontent.com/u/15937715/Data/github/goatools/gograph-gml.png
    :alt: GO term lineage (Cytoscape)


Map GO terms to GOslim terms
::::::::::::::::::::::::::::::::::::
See ``map_to_slim.py`` for usage. As arguments it takes the gene ontology files:

* the current gene ontology file ``go-basic.obo``
* the GOslim file to be used (e.g. ``goslim_generic.obo`` or any other GOslim
  file)

The script either maps one GO term to it's GOslim terms, or protein products
with multiple associations to all it's GOslim terms.

To determine the GOslim terms for a single GO term, you can use the following
command::

    python scripts/map_to_slim.py --term=GO:0008135 go-basic.obo goslim_generic.obo

To determine the GOslim terms for protein products with multiple associations::

    python scripts/map_to_slim.py --association_file=data/association go-basic.obo goslim_generic.obo

Where the ``association`` file has the same format as used for
``find_enrichment.py``.

The implemented algorithm is described in more detail at the go-perl
documenation of `map2slim <http://search.cpan.org/~cmungall/go-perl/scripts/map2slim>`_.
