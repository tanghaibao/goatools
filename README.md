# Tools for Gene Ontology

[![DOI](/doc/images/DOI.svg)](https://www.nature.com/articles/s41598-018-28948-z)
[![Latest PyPI version](https://img.shields.io/pypi/v/goatools.svg)](https://pypi.python.org/pypi/goatools)
[![bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/goatools/README.html?highlight=goatools)
[![Github Actions](https://github.com/tanghaibao/goatools/workflows/build/badge.svg)](https://github.com/tanghaibao/goatools/actions)
[![Downloads](https://pepy.tech/badge/goatools)](https://pepy.tech/project/goatools)

|         |                                                                       |
| ------- | --------------------------------------------------------------------- |
| Authors | Haibao Tang ([tanghaibao](http://github.com/tanghaibao))              |
|         | DV Klopfenstein ([dvklopfenstein](https://github.com/dvklopfenstein)) |
|         | Brent Pedersen ([brentp](http://github.com/brentp))                   |
|         | Fidel Ramirez ([fidelram](https://github.com/fidelram))               |
|         | Aurelien Naldi ([aurelien-naldi](http://github.com/aurelien-naldi))   |
|         | Patrick Flick ([patflick](http://github.com/patflick))                |
|         | Jeff Yunes ([yunesj](http://github.com/yunesj))                       |
|         | Kenta Sato ([bicycle1885](http://github.com/bicycle1885))             |
|         | Chris Mungall ([cmungall](https://github.com/cmungall))               |
|         | Greg Stupp ([stuppie](https://github.com/stuppie))                    |
|         | David DeTomaso ([deto](https://github.com/deto))                      |
|         | Olga Botvinnik ([olgabot](https://github.com/olgabot))                |
| Email   | <tanghaibao@gmail.com>                                                |
| License | BSD                                                                   |

## Description

This package contains a Python library to

- Process over- and under-representation of certain GO terms, based on
  Fisher's exact test. With numerous multiple correction routines
  including locally implemented routines for Bonferroni, Sidak, Holm,
  and false discovery rate. Also included are multiple test
  corrections from
  [statsmodels](http://www.statsmodels.org/stable/index.html): FDR
  Benjamini/Hochberg, FDR Benjamini/Yekutieli, Holm-Sidak,
  Simes-Hochberg, Hommel, FDR 2-stage Benjamini-Hochberg, FDR 2-stage
  Benjamini-Krieger-Yekutieli, FDR adaptive Gavrilov-Benjamini-Sarkar,
  Bonferroni, Sidak, and Holm.

- Process the obo-formatted file from [Gene Ontology
  website](http://geneontology.org). The data structure is a directed
  acyclic graph (DAG) that allows easy traversal from leaf to root.

- Read [GO Association files](http://geneontology.org/page/go-annotation-file-formats):

  - GAF ([GO Annotation
    File](http://geneontology.org/page/go-annotation-file-gaf-format-21))
  - GPAD ([Gene Product Association Data](https://geneontology.github.io/docs/gene-product-association-data-gpad-format/))
  - NCBI's gene2go file
  - id2gos format. See [example](https://raw.githubusercontent.com/tanghaibao/goatools/main/data/association)

- [Print **_decendants count_** and/or **_information content_**](/notebooks/dcnt_and_tinfo.ipynb)
  for a list of GO terms

- [Get parents or ancestors for a GO term with or without optional
  relationships](notebooks/parents_and_ancestors.ipynb), including [Print details about a
  GO ID's parents](notebooks/parent_go_terms.ipynb)

- Compare two or more lists of GO IDs
- [Plot GO hierarchies](https://github.com/tanghaibao/goatools#plot-go-lineage)
- [Write GO hierarchies to an ASCII text file](https://github.com/tanghaibao/goatools#write-go-hierarchy)
- Group GO terms for easier viewing

- Map GO terms (or protein products with multiple associations to
  GO terms) to GOslim terms (analog to the map2slim.pl script supplied
  by geneontology.org)

## To Cite

Klopfenstein DV, Zhang L, Pedersen BS, ... Tang H [GOATOOLS: A Python library for Gene Ontology analyses](https://www.nature.com/articles/s41598-018-28948-z)
_Scientific reports_ | (2018) 8:10872 | DOI:10.1038/s41598-018-28948-z

- **GO Grouping**:
  Visualize the major findings in a gene ontology enrichment analysis (GEOA) more easily with grouping.
  A detailed description of GOATOOLS GO grouping is found in the
  [manuscript](https://www.nature.com/articles/s41598-018-28948-z).
- **Compare GO lists**:
  Compare [two](https://github.com/tanghaibao/goatools/issues/162) or more lists
  of GO IDs using `compare_gos.py`, which can be used with or without
  grouping.
- **Stochastic GOEA simulations**:
  One of the findings resulting from our simulations is:
  [Larger study sizes result in higher GOEA sensitivity](https://github.com/dvklopfenstein/goatools_simulation#manuscript-figures),
  meaning fewer truly significant observations go unreported.
  The code for the stochastic GOEA simulations
  described in the paper is found [here](https://github.com/dvklopfenstein/goatools_simulation)

## Installation

Make sure your Python version >= 3.7, install the latest stable
version via PyPI:

```bash
pip install goatools
```

To install the development version:

```bash
pip install git+git://github.com/tanghaibao/goatools.git
```

`.obo` file for the most current
[GO](http://geneontology.org/page/download-ontology):

```bash
wget http://geneontology.org/ontology/go-basic.obo
```

`.obo` file for the most current [GO
Slim](http://geneontology.org/page/go-slim-and-subset-guide) terms (e.g.
generic GOslim) :

```bash
wget http://www.geneontology.org/ontology/subsets/goslim_generic.obo
```

## Dependencies

- Simplest is to install via pip as above or
  [bioconda](http://bioconda.github.io/recipes/goatools/README.html?highlight=goatools).

- For statistical testing of GO enrichment:

  - [`scipy.stats.fisher_exact`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html)
  - `statsmodels` (optional) for access to a variety of statistical tests for GOEA

- To plot the ontology lineage, install one of these two options:

  - [Graphviz](http://www.graphviz.org/), for graph visualization.
  - [pygraphviz](http://networkx.lanl.gov/pygraphviz/), Python binding for communicating with Graphviz:
  - [pydot](https://code.google.com/p/pydot/), a Python interface to Graphviz's Dot language.

## Cookbook

`run.sh` contains example cases, which calls the utility scripts in the
`scripts` folder.

### Find GO enrichment of genes under study

See examples in [find_enrichment](/doc/md/README_find_enrichment.md)

The `find_enrichment.py` takes as arguments files
containing:

- gene names in a study
- gene names in population (or other study if `--compare` is specified)
- an association file that maps a gene name to a GO category.

Please look at `tests/data` folder to see examples on how to make these
files. when ready, the command looks like:

```bash
python scripts/find_enrichment.py --pval=0.05 --indent data/study \
                                  data/population data/association
```

and can filter on the significance of (e)nrichment or (p)urification. it
can report various multiple testing corrected p-values as well as the
false discovery rate.

The `e` in the "Enrichment" column means "enriched" - the concentration
of GO term in the study group is significantly _higher_ than those in
the population. The "p" stands for "purified" - significantly _lower_
concentration of the GO term in the study group than in the population.

**Important note**: by default, `find_enrichment.py` propagates counts
to all the parents of a GO term. As a result, users may find terms in
the output that are not present in their `association` file. Use
`--no_propagate_counts` to disable this behavior.

### Write GO hierarchy

- [`wr_hier.py`](doc/md/README_wr_hier.md): Given a GO ID, write the hierarchy below (default)
  or above (`--up`) the given GO.

### Plot GO lineage

- [`go_plot.py`](doc/md/README_go_plot.md):
  - Plots user-specified GO term(s) up to root
  - Multiple user-specified GOs
  - User-defined colors
  - Plot relationships (`-r`)
  - Optionally plot children of user-specfied GO terms
- `plot_go_term.py` can plot the lineage
of a certain GO term, by:

```bash
python scripts/plot_go_term.py --term=GO:0008135
```

This command will plot the following image.

![GO term lineage](https://www.dropbox.com/s/4zbqx8sqcls3mge/gograph.png?raw=1)

Sometimes people like to stylize the graph themselves, use option
`--gml` to generate a GML output which can then be used in an external
graph editing software like [Cytoscape](http://www.cytoscape.org/). The
following image is produced by importing the GML file into Cytoscape
using yFile orthogonal layout and solid VizMapping. Note that the [GML
reader plugin](https://code.google.com/p/graphmlreader/) may need to be
downloaded and installed in the `plugins` folder of Cytoscape:

```bash
python scripts/plot_go_term.py --term=GO:0008135 --gml
```

![GO term lineage (Cytoscape)](https://www.dropbox.com/s/ueov2ioxl063q8h/gograph-gml.png?raw=1)

### Map GO terms to GOslim terms

See `map_to_slim.py` for usage. As arguments it takes the gene ontology
files:

- the current gene ontology file `go-basic.obo`
- the GOslim file to be used (e.g. `goslim_generic.obo` or any other GOslim file)

The script either maps one GO term to its GOslim terms, or protein
products with multiple associations to all its GOslim terms.

To determine the GOslim terms for a single GO term, you can use the
following command:

```bash
python scripts/map_to_slim.py --term=GO:0008135 go-basic.obo goslim_generic.obo
```

To determine the GOslim terms for protein products with multiple
associations:

```bash
python scripts/map_to_slim.py --association_file=data/association go-basic.obo goslim_generic.obo
```

Where the `association` file has the same format as used for
`find_enrichment.py`.

The implementation is similar to
[map2slim](http://search.cpan.org/~cmungall/go-perl/scripts/map2slim).

## Technical notes

### Available statistical tests for calculating uncorrected p-values

For calculating uncorrected p-values, we use SciPy:

- [`scipy.stats.fisher_exact`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html)

### Available multiple test corrections

We have implemented several significance tests:

- `bonferroni`, bonferroni correction
- `sidak`, sidak correction
- `holm`, hold correction
- `fdr`, false discovery rate (fdr) implementation using resampling

Additional methods are available if `statsmodels` is installed:

- `sm_bonferroni`, bonferroni one-step correction
- `sm_sidak`, sidak one-step correction
- `sm_holm-sidak`, holm-sidak step-down method using Sidak adjustments
- `sm_holm`, holm step-down method using Bonferroni adjustments
- `simes-hochberg`, simes-hochberg step-up method (independent)
- `hommel`, hommel closed method based on Simes tests (non-negative)
- `fdr_bh`, fdr correction with Benjamini/Hochberg (non-negative)
- `fdr_by`, fdr correction with Benjamini/Yekutieli (negative)
- `fdr_tsbh`, two stage fdr correction (non-negative)
- `fdr_tsbky`, two stage fdr correction (non-negative)
- `fdr_gbs`, fdr adaptive Gavrilov-Benjamini-Sarkar

In total 15 tests are available, which can be selected using option
`--method`. Please note that the default FDR (`fdr`) uses a resampling
strategy which may lead to slightly different q-values between runs.

## iPython [Notebooks](notebooks)

### Optional attributes

[definition](notebooks/Optional_attribute_definition.ipynb)

### Run a Ontology Enrichment Analysis (GOEA)

[goea_nbt3102](notebooks/goea_nbt3102.ipynb)    
[human phenotype ontologies](notebooks/Enrichment_analyses_human_phenotype_ontology.ipynb)    

### Show many study genes are associated with RNA, translation, mitochondria, and ribosomal

[goea_nbt3102_group_results](notebooks/goea_nbt3102_group_results.ipynb)

### Report level and depth counts of a set of GO terms

[report_depth_level](notebooks/report_depth_level.ipynb)

### Find all human protein-coding genes associated with cell cycle

[cell_cycle](notebooks/cell_cycle.ipynb)

### Calculate annotation coverage of GO terms on various species

[annotation_coverage](notebooks/annotation_coverage.ipynb)

### Determine the semantic similarities between GO terms

[semantic_similarity](notebooks/semantic_similarity.ipynb) [semantic_similarity_wang](notebooks/semantic_similarity_wang.ipynb)

### Obsolete GO terms are loaded upon request

[godag_obsolete_terms](notebooks/godag_obsolete_terms.ipynb)

## Want to Help?

Prior to submitting your pull request, please add a test which verifies your code, and run:

```console
make test
```

Items that we know we need include:

- Add code coverage runs
- Edit tests in the `makefile` under the comment
- Help setting up [documentation](http://goatools.readthedocs.io/en/latest/). We
  are using Sphinx and Python docstrings to create documentation.
  For documentation practice, use make targets:

  ```bash
  make mkdocs_practice
  ```

  To remove practice documentation:

  ```bash
  make rmdocs_practice
  ```

  Once you are happy with the documentation do:

  ```bash
  make gh-pages
  ```

Copyright (C) 2010-2021, Haibao Tang et al. All rights reserved.
