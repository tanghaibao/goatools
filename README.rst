Tools for Gene Ontology
========================

.. image:: https://zenodo.org/badge/doi/10.5281/zenodo.31628.svg
    :target: http://dx.doi.org/10.5281/zenodo.31628
    :alt: DOI

.. image:: https://img.shields.io/pypi/v/goatools.svg
    :target: https://crate.io/packages/goatools/
    :alt: Latest PyPI version

.. image:: https://img.shields.io/pypi/dm/goatools.svg
    :target: https://crate.io/packages/goatools/
    :alt: Number of PyPI downloads

.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat
    :target: http://bioconda.github.io/recipes/goatools/README.html?highlight=goatools
    :alt: bioconda

.. image:: https://travis-ci.org/tanghaibao/goatools.svg?branch=master
    :target: https://travis-ci.org/tanghaibao/goatools
    :alt: Travis-CI

:Author: Haibao Tang (`tanghaibao <http://github.com/tanghaibao>`_),
         Brent Pedersen (`brentp <http://github.com/brentp>`_),
         Fidel Ramirez (`fidelram <https://github.com/fidelram>`_),
         Aurelien Naldi (`aurelien-naldi <http://github.com/aurelien-naldi>`_),
         Patrick Flick (`patflick <http://github.com/patflick>`_),
         Jeff Yunes (`yunesj <http://github.com/yunesj>`_),
         Kenta Sato (`bicycle1885 <http://github.com/bicycle1885>`_),
         Chris Mungall (`cmungall <https://github.com/cmungall>`_),
         Greg Stupp (`stuppie <https://github.com/stuppie>`_),
         DV Klopfenstein(`dvklopfenstein <https://github.com/dvklopfenstein>`_)
:Email: tanghaibao@gmail.com
:License: BSD

.. contents ::

Description
------------
This package contains a Python library to

- process over- and under-representation of certain GO terms, based on Fisher's
  exact test. With numerous multiple correction routines including locally
  implemented routines for Bonferroni, Sidak, Holm, and false discovery rate. Also included are
  multiple test corrections from `statsmodels <http://www.statsmodels.org/stable/index.html>`_:
  FDR Benjamini/Hochberg, FDR Benjamini/Yekutieli, Holm-Sidak, Simes-Hochberg,
  Hommel, FDR 2-stage Benjamini-Hochberg, FDR 2-stage Benjamini-Krieger-Yekutieli,
  FDR adaptive Gavrilov-Benjamini-Sarkar, Bonferroni, Sidak, and Holm.
- process the obo-formatted file from `Gene Ontology website <http://geneontology.org>`_.
  The data structure is a directed acyclic graph (DAG) that allows easy traversal
  from leaf to root.
- map GO terms (or protein products with multiple associations to GO terms) to
  GOslim terms (analog to the map2slim.pl script supplied by geneontology.org)


Installation
-------------
Make sure your Python version >= 2.7, install the latest stable version via PyPI::

    easy_install goatools

To install the development version::

    pip install git+git://github.com/tanghaibao/goatools.git

``.obo`` file for the most current `GO <http://geneontology.org/page/download-ontology>`_::

    wget http://geneontology.org/ontology/go-basic.obo

``.obo`` file for the most current `GO Slim <http://geneontology.org/page/go-slim-and-subset-guide>`_
terms (e.g. generic GOslim) ::

    wget http://www.geneontology.org/ontology/subsets/goslim_generic.obo


Dependencies
-------------
- `fisher <http://pypi.python.org/pypi/fisher/>`_ (required) for calculating
  Fisher's exact test::

    easy_install fisher

- `statsmodels` (optional) for access to a variety of statistical tests for GOEA::

    easy_install statsmodels

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

- Alternatively, it is possible to install via bioconda. See details
`here <http://bioconda.github.io/recipes/goatools/README.html?highlight=goatools>`_.


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

**Important note**: by default, ``find_enrichment.py`` propagates counts to all
the parents of a GO term. As a result, users may find terms in the output that
are not present in their ``association`` file. Use ``--no_propagate_counts`` to
disable this behavior.

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

The script either maps one GO term to its GOslim terms, or protein products
with multiple associations to all its GOslim terms.

To determine the GOslim terms for a single GO term, you can use the following
command::

    python scripts/map_to_slim.py --term=GO:0008135 go-basic.obo goslim_generic.obo

To determine the GOslim terms for protein products with multiple associations::

    python scripts/map_to_slim.py --association_file=data/association go-basic.obo goslim_generic.obo

Where the ``association`` file has the same format as used for
``find_enrichment.py``.

The implemented algorithm is described in more detail at the go-perl
documentation of `map2slim <http://search.cpan.org/~cmungall/go-perl/scripts/map2slim>`_.


Available significance tests
::::::::::::::::::::::::::::
We have implemented several significance tests:

* ``bonferroni``, bonferroni correction
* ``sidak``, sidak correction
* ``holm``, hold correction
* ``fdr``, false discovery rate (fdr) implementation using resampling

Additional methods are available if ``statsmodels`` is installed:

* ``sm_bonferroni``, bonferroni one-step correction
* ``sm_sidak``, sidak one-step correction
* ``sm_holm-sidak``, holm-sidak step-down method using Sidak adjustments
* ``sm_holm``, holm step-down method using Bonferroni adjustments
* ``simes-hochberg``, simes-hochberg step-up method (independent)
* ``hommel``, hommel closed method based on Simes tests (non-negative)
* ``fdr_bh``, fdr correction with Benjamini/Hochberg (non-negative)
* ``fdr_by``, fdr correction with Benjamini/Yekutieli (negative)
* ``fdr_tsbh``, two stage fdr correction (non-negative)
* ``fdr_tsbky``, two stage fdr correction (non-negative)
* ``fdr_gbs``, fdr adaptive Gavrilov-Benjamini-Sarkar

In total 15 tests are available, which can be selected using option ``--method``.
Please note that the default FDR (``fdr``) uses a resampling strategy which may
lead to slightly different q-values between runs.


iPython Notebooks
-----------------

Run a Gene Ontology Enrichment Analysis (GOEA)
::::::::::::::::::::::::::::::::::::::::::::::
https://github.com/tanghaibao/goatools/blob/master/notebooks/goea_nbt3102.ipynb

Show many study genes are associated with RNA, translation, mitochondria, and ribosomal
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
https://github.com/tanghaibao/goatools/blob/master/notebooks/goea_nbt3102_group_results.ipynb

Report level and depth counts of a set of GO terms
::::::::::::::::::::::::::::::::::::::::::::::::::
https://github.com/tanghaibao/goatools/blob/master/notebooks/report_depth_level.ipynb

Find all human protein-coding genes associated with cell cycle
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
https://github.com/tanghaibao/goatools/blob/master/notebooks/cell_cycle.ipynb

Calculate annotation coverage of GO terms on various species
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
https://github.com/tanghaibao/goatools/blob/master/notebooks/annotation_coverage.ipynb

Determine the semantic similarities between GO terms
::::::::::::::::::::::::::::::::::::::::::::::::::::
https://github.com/tanghaibao/goatools/blob/master/notebooks/semantic_similarity.ipynb


Reference
---------
Haibao Tang et al. (2015). GOATOOLS: Tools for Gene Ontology. Zenodo.
`10.5281/zenodo.31628 <http://dx.doi.org/10.5281/zenodo.31628>`_.
