# CHANGELOG

## Summary

* [**Unreleased changes**](#unreleased-changes)
* [**Release 2021-05-22 1.1.5**](#release-2021-05-22-115)
* [**Release 2020-12-04 1.0.13**](#release-2020-12-04-1013)
* [**Release 2020-06-23 1.0.3**](#release-2020-06-23-106)
  * Added code to download background genes from NCBI.
  * Fixed the 'write hierarchy' script, wr_hier.py, so it will write MF and CC hierarchies.
* [**Release 2020-03-13 1.0.3**](#release-2020-03-13-103) GOEAs can be run in quiet mode
* [**Release 2020-02-20 1.0.2**](#release-2020-02-20-102)
  * Added [Jupyter notebook](https://github.com/tanghaibao/goatools/blob/main/notebooks/godag_obsolete_terms.ipynb) showing how to work with obsolete GO terms
    [#153](https://github.com/tanghaibao/goatools/issues/153)
    [#154](https://github.com/tanghaibao/goatools/issues/154)
    [#155](https://github.com/tanghaibao/goatools/issues/155)
  * Deprecated: Internal data member, *go2parents* will be deprecated, renamed to *go2ancestors*
* [**Release 2019-09-29 0.9.7**](#release-2019-09-29-097)
  * Deprecated: *read_ncbi_gene2go*
  * Deprecated: *get_b2aset* and *GoDagTimed* in their old location. They have been moved.
  * Removed empty sets from Descendant and Ancestor dicts in *gosubdag.rcntobj*
  * *TermCount* improvements (used in semantic similarity calculations)
  * GO DAG Plotting:
    * User can now add a title to GO DAG plots
    * User can now add custom text to edges in GO DAG plots
* [**Release 2019-07-26 0.9.7**](#release-2019-07-26-097)
  * Support user-specified optional relationships in GOEA propagate counts
* [**Release 2019-05-08 0.9.5**](#release-2019-05-08-095)
  * Support user-specified evidence codes in GOEAs
  * Separate GOEAs into BP, MF, and CC

## Unreleased changes
* **Added**
  * Added human phenotype ontologies enrichement analyses [#202](https://github.com/tanghaibao/goatools/issues/202)

## Release 2021-05-22 1.1.5

* **Added**
  * Support for GAF 2.2 [#201](https://github.com/tanghaibao/goatools/issues/201)
* **Changed**
  * The "Compare GO terms" script now can use GO terms stored in files of any name. [#187](https://github.com/tanghaibao/goatools/issues/187)
  * Improved documentation about using the optional attribute, definition. [#204](https://github.com/tanghaibao/goatools/issues/204)

## Release 2020-12-04 1.0.13

* **Added**
  * Added Wang's semantic similarity calculator
  * Added a script to duplicate Table 2 in the [GOATOOLS publication](https://www.nature.com/articles/s41598-018-28948-z) [#171](https://github.com/tanghaibao/goatools/issues/171)
  * Added clear documentation to explain GO term text (D and R) in plots [#178](https://github.com/tanghaibao/goatools/issues/178)
* **Changed**
  * Workarounds for errors in the Gene Ontology Consortium (GOC) annotations and opened these issues with GOC.
    NOTE: Report annotation issues to [Help Desk](https://github.com/geneontology/helpdesk/)
    * [Incorrect taxon format](https://github.com/geneontology/go-annotation/issues/3373) [GO helpdesk #287](https://github.com/geneontology/helpdesk/issues/287)
    * [Obsolete GOs seen in gaf files](https://github.com/geneontology/go-annotation/issues/3523)
  * Moved calculating the depth using optional relationships from GODag to its own module [#188](https://github.com/tanghaibao/goatools/issues/188)

## Release 2020-06-23 1.0.6

* **Fixed**
  * The `scripts/wr_hier.pyi` script can now print the hierarchy for GO IDs in all namespaces. [#163](https://github.com/tanghaibao/goatools/issues/163)
  * The semantic_similarity function upon comparing a GO term to itself. [#150](https://github.com/tanghaibao/goatools/issues/150)
* **Added**
  * Added a new notebook showing how to download background genes from NCBI.
  * Added arg, `--prt_study_gos_only`, to script, `scripts/find_enrichment.py`
    to print only study GOs when printing all GO terms, regardless of their significance (`--pval=1.0`):
    `find_enrichment.py study_genes.txt human_genes.txt gene2go --pval=1.0 --prt_study_gos_only`
* **Changed**
  * Remove trailing divider ("NOT|"), if it exists in the gpad file ([go-annotation #2885](https://github.com/geneontology/go-annotation/issues/2885))

## Release 2020-03-13 1.0.3

* **Added**
  * Added quiet mode when running GOEAs using the function, *run_goea*.
    Examples are in section **5a. Quietly Run a GOEA** of the
[**GOEA Jupyter notebook**](https://github.com/tanghaibao/goatools/blob/main/notebooks/goea_nbt3102.ipynb).
re: See comments ([1](https://github.com/tanghaibao/goatools/issues/133#issuecomment-598334391) and
[2](https://github.com/tanghaibao/goatools/issues/133#issuecomment-598844019)) in
[**#133**](https://github.com/tanghaibao/goatools/issues/133)

## Release 2020-02-20 1.0.2

* **Added**
  * New [Jupyter notebook](https://github.com/tanghaibao/goatools/blob/main/notebooks/godag_obsolete_terms.ipynb) showing how to work with obsolete GO terms.
    [**#153**](https://github.com/tanghaibao/goatools/issues/153)
    [**#154**](https://github.com/tanghaibao/goatools/issues/154)
    [**#155**](https://github.com/tanghaibao/goatools/issues/155)
  * [**Issue 154**](https://github.com/tanghaibao/goatools/issues/154) Support for optional GO term attributes, *consider* and *replaced_by*
* **Deprecated**
  * Renamed internal data variable:
    * NOW: *gosubdag.rcntobj.go2ancestors*
    * WAS: *gosubdag.rcntobj.go2parents*
* **Changed**
  * [**Issue 148**](https://github.com/tanghaibao/goatools/issues/148)
    Return None for Lin's semantic similarity calculations if one or both of the GO terms is not annotated.
    This will be updated in the future to have the option to assign a count of 1 to GO terms that are not annotated,
    indicating the annotation of a mock gene to allow the researcher to get a rough idea of the similarity.
  * [**Issue 142**](https://github.com/tanghaibao/goatools/issues/142)
    * Write GO hierachy to a file now writes a file when using Python3
    * Gaf reader defaults to gaf file version of 2.1 if no version line if found

## Release 2019-09-29 0.9.9

* **Deprecated**
  * *read_ncbi_gene2go* is deprecated and will be removed in the future.
  * *get_b2aset* is moved:
    * NOW: goatools.utils
    * WAS: goatools.associations
  * *GoDagTimed* and *prt_hms* is moved:
    * NOW: goatools.godag.prttime
    * WAS: goatools.test_data.godag_timed
* **Removed**
  * All dict entries whose values were an empty set in:
    * *gosubdag.rcntobj.go2parents*
    * *gosubdag.rcntobj.go2descendants*
* **Added**
  * Method to annotation object, *IdToGosReader*, which writes namedtuples into an ASCII file
* **Changed**
  * *TermCounts*:
    * Added support for optional relationships, like *part_of*.
      This is useful for computing termwise and genewise semantic similarities.
    * Faster initialization of *TermCounts* object, used in semantic similarity calculations
  * *Plotting*:
    * Users can now provide a title to be printed in a GO DAG plot
    * Users can now provide an *edge2txt* dict to print text on edges between GO Terms
* **Fixed**
  * Aspect counts (BP, MF, CC totals) in *TermCounts* object explained in [#156](https://github.com/tanghaibao/goatools/issues/156)
  * Opened issue:
    [#2629](https://github.com/geneontology/go-annotation/issues/2629) -> [Help
    Desk](https://github.com/geneontology/helpdesk/issues/288)

## Release 2019-07-26 0.9.7

* **Added**
  * Support traversing optional relationships when propagating counts for GOEAs
  * Find all ancesters of a GO term using a user-specified list of relationships
  [#126](https://github.com/tanghaibao/goatools/issues/126#issuecomment-511985524)

### Support traversing optional relationships

Support traversing optional relationships, like *part_of* and *regulates*,
when doing propagating counts using these two new arguments in `scripts/find_enrichment.py`:

```bash
        -r, --relationship    Propagate counts up all relationships (default: False)
       --relationships        [RELATIONSHIPS [RELATIONSHIPS ...]]
                               Propagate counts up user-specified relationships
                               (default: None)
```

[#117](https://github.com/tanghaibao/goatools/issues/117#issuecomment-515484624)

## Release 2019-05-08 0.9.5

* **Added**
  * [**Issue 119: Added support for specifying specific evidence codes**](#added-support-for-specifying-specific-evidence-codes)
  * [**Issue 127: Added splitting GOEA into BP, MF, and CC**](#added-splitting-goea-into-bp-mf-and-cc)
* **Fixed**
  * [**Issue 120 and 121: Lin similarity is positive**](https://github.com/tanghaibao/goatools/issues/120)

### Added support for specifying specific evidence codes

* Specify evidence codes, like EXP (Inferred from Experiment), to exclude or include in a GOEA.
* Specify evidence classes, like Experimental (EXP IDA IPI IMP IGI IEP), which include many evidence codes.
* Get evidence code help:

```bash
python3 scripts/find_enrichment.py --ev_help
python3 scripts/find_enrichment.py --ev_help_short
```

[#119](https://github.com/tanghaibao/goatools/issues/119)

### Added splitting GOEA into BP, MF, and CC

Split GOEA into three separate analyses by default:

* BP (biological process)
* MF (molecular function)
* CC (cellular component)

[#127](https://github.com/tanghaibao/goatools/issues/127#issuecomment-489776548)

## Release 2018-11-27 0.8.12

## Release 2018-11-21 0.8.11

## Release 2018-09-09 0.8.9

## Release 2018-04-28 0.8.4

## Release 2018-02-22 0.8.2

## Release 2018-02-19 0.7.11

## Release 2018-01-07 0.7.11

## Release 2017-10-14 0.7.11

## Release 2017-10-07 0.7.10

## Release 2015-02-05

## Release 2014-08-10

## Release 2010-04-12

## Release 2010-02-23

[keep a changelog](https://keepachangelog.com/en/1.0.0/)
