# Tools for Gene Ontology

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.31628.svg)](http://dx.doi.org/10.5281/zenodo.31628)
[![Latest PyPI version](https://img.shields.io/pypi/v/goatools.svg)](https://pypi.python.org/pypi/goatools)
[![Number of PyPI downloads](https://img.shields.io/pypi/dm/goatools.svg)](https://pypi.python.org/pypi/goatools)
[![bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/goatools/README.html?highlight=goatools)
[![Travis-CI](https://travis-ci.org/tanghaibao/goatools.svg?branch=master)](https://travis-ci.org/tanghaibao/goatools)

| | |

90909090909
llllllllllllllllllllllllllll


ooooooooooooooooooooooooooooooooooollllllllllllllllllllllllllllllllllllllllllll
|---|---|
| Author | Haibao Tang ([tanghaibao](http://github.com/tanghaibao)) |
| | DV Klopfenstein ([dvklopfenstein](https://github.com/dvklopfenstein)) |
| | Brent Pedersen ([brentp](http://github.com/brentp)) |
| | Fidel Ramirez ([fidelram](https://github.com/fidelram)) |
| | Aurelien Naldi ([aurelien-naldi](http://github.com/aurelien-naldi)) |
| | Patrick Flick ([patflick](http://github.com/patflick)) |
| | Jeff Yunes ([yunesj](http://github.com/yunesj)) |
| | Kenta Sato ([bicycle1885](http://github.com/bicycle1885)) |
| | Chris Mungall ([cmungall](https://github.com/cmungall)) |
| | Greg Stupp ([stuppie](https://github.com/stuppie)) |
| | David DeTomaso ([deto](https://github.com/deto)) |
| | Olga Botvinnik ([olgabot](https://github.com/olgabot)) |
| Email | <tanghaibao@gmail.com> |
| License | BSD |

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

- Read [GO Association
  files](http://geneontology.org/page/go-annotation-file-formats):
  - Read GAF ([Gene Association
      File](http://geneontology.org/page/go-annotation-file-gaf-format-21)) files.
  - Read NCBI's gene2go GO association file.

- Map GO terms (or protein products with multiple associations to
  GO terms) to GOslim terms (analog to the map2slim.pl script supplied
  by geneontology.org)

## Installation

Make sure your Python version >= 2.7, install the latest stable
version via PyPI:

```bash
easy_install goatools
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

- Simplest is to install via bioconda. See details
   [here](http://bioconda.github.io/recipes/goatools/README.html?highlight=goatools).

- To calculate the uncorrected p-values, there are currently twooptions:
  - [fisher](http://pypi.python.org/pypi/fisher/) for calculating Fisher's exact test:

  ```bash
  easy_install fisher
  ```

  - [fisher](http://docs.scipy.org/doc/scipy-0.17.0/reference/generated/scipy.statfisher_exact.html)
    from [SciPy's](http://docs.scipy.org/doc/scipy/reference/)
    [stats](http://docs.scipy.org/doc/scipy/reference/stats.html) package

  - `statsmodels` (optional) for access to a variety of statistical tests for GOEA:

   ```bash
   easy_install statsmodels
   ```

- To plot the ontology lineage, install one of these two options:
  - Graphviz
    - [Graphviz](http://www.graphviz.org/), for graph visualization.
    - [pygraphviz](http://networkx.lanl.gov/pygraphviz/), Python binding for communicating with Graphviz:

    ```bash
    easy_install pygraphviz
    ```

  - [pydot](https://code.google.com/p/pydot/), a Python interface to Graphviz's Dot language.
    - [pyparsing](http://pyparsing.wikispaces.com/) is a prerequisite for `pydot`
    - Images can be viewed using either:
      - [ImageMagick](http://www.imagemagick.org/)'s *display*
      - [Graphviz](http://www.graphviz.org/)

## Cookbook

`run.sh` contains example cases, which calls the utility scripts in the
`scripts` folder.

### Find GO enrichment of genes under study

See `find_enrichment.py` for usage. It takes as arguments files
containing:

- gene names in a study
- gene names in population (or other study if `--compare` is specified)
- an association file that maps a gene name to a GO category.

Please look at `tests/data/` folder to see examples on how to make these
files. when ready, the command looks like:

```bash
python scripts/find_enrichment.py --pval=0.05 --indent data/study \
                                  data/population data/association
```

and can filter on the significance of (e)nrichment or (p)urification. it
can report various multiple testing corrected p-values as well as the
false discovery rate.

The "e" in the "Enrichment" column means "enriched" - the concentration
of GO term in the study group is significantly *higher* than those in
the population. The "p" stands for "purified" - significantly *lower*
concentration of the GO term in the study group than in the population.

**Important note**: by default, `find_enrichment.py` propagates counts
to all the parents of a GO term. As a result, users may find terms in
the output that are not present in their `association` file. Use
`--no_propagate_counts` to disable this behavior.

### Read and plot GO lineage

See `plot_go_term.py` for usage. `plot_go_term.py` can plot the lineage
of a certain GO term, by:

```bash
python scripts/plot_go_term.py --term=GO:0008135
```

This command will plot the following image.

![GO term lineage](https://dl.dropboxusercontent.com/u/15937715/Data/github/goatools/gograph.png)

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

![GO term lineage (Cytoscape)](https://dl.dropboxusercontent.com/u/15937715/Data/github/goatools/gograph-gml.png)

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

The implemented algorithm is described in more detail at the go-perl
documentation of
[map2slim](http://search.cpan.org/~cmungall/go-perl/scripts/map2slim).

## Technical notes

### Available statistical tests for calculating uncorrected p-values

There are currently two fisher tests available for calculating uncorrected
p-values. Both fisher options from the fisher package and SciPy's stats package
calculate the same pvalues, but provide the user an option in installing
packages.

- `fisher`, [fisher](http://pypi.python.org/pypi/fisher/) package's `fisher.pvalue_population`
- `fisher_scipy_stats`:[SciPy](http://docs.scipy.org/doc/scipy/reference)
   [stats](http://docs.scipy.org/doc/scipy/reference/stats.html) package
  [fisher_exact](http://docs.scipy.org/doc/scipy-0.17.0/reference/generated/scipy.statfisher_exact.html)

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

## iPython Notebooks

### Run a Gene Ontology Enrichment Analysis (GOEA)

<https://github.com/tanghaibao/goatools/blob/master/notebooks/goea_nbt3102.ipynb>

### Show many study genes are associated with RNA, translation, mitochondria, and ribosomal

<https://github.com/tanghaibao/goatools/blob/master/notebooks/goea_nbt3102_group_results.ipynb>

### Report level and depth counts of a set of GO terms

<https://github.com/tanghaibao/goatools/blob/master/notebooks/report_depth_level.ipynb>

### Find all human protein-coding genes associated with cell cycle

<https://github.com/tanghaibao/goatools/blob/master/notebooks/cell_cycle.ipynb>

### Calculate annotation coverage of GO terms on various species

<https://github.com/tanghaibao/goatools/blob/master/notebooks/annotation_coverage.ipynb>

### Determine the semantic similarities between GO terms

<https://github.com/tanghaibao/goatools/blob/master/notebooks/semantic_similarity.ipynb>

## Want to Help?

If you add new code, please be sure to also add python tests which
verify your code.

Items that we know we need include:

- Add code coverage runs
- Edit tests in the `makefile` under the comment, `# TBD`, suchthey run using `nosetests`
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

## Reference

Haibao Tang et al. (2015). GOATOOLS: Tools for Gene Ontology. Zenodo.
[10.5281/zenodo.31628](http://dx.doi.org/10.5281/zenodo.31628).
