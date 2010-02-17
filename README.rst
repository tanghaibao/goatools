Description
-----------
This package contains two utility scripts to handle gene ontology enrichment analysis.


Installation
------------
The only requirement is to run the `setup.py` script::

    python setup.py build_ext -i

this builds the C-extension of Fisher's exact test for speed

Scripts
-------
``genemerge.py``
    command-line tool for processing over- and under-representation of certain GO terms, based on Fisher's exact test. Implemented several multiple correction routines (including Bonferroni, sidak, and false discovery rate).

``obo_parser.py``
    simple library to process the obo-formatted file from `Gene Ontology website <http://geneontology.org>`_. The data structure is a directed acyclic graph (DAG) that allows easy traversal from leaf to root.

