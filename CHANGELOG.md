# CHANGELOG

### Summary

* [**Unreleased changes**](#unreleased-changes)
  * Deprecated: Internal data member, *go2parents* will be deprecated, renamed to *go2ancestors*
  * Changes for [**Issue 142**](https://github.com/tanghaibao/goatools/issues/142)    
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

Unreleased changes
--------------

* **Deprecated**
  * Renamed internal data variable:
    * NOW: *gosubdag.rcntobj.go2ancestors*    
    * WAS: *gosubdag.rcntobj.go2parents*    
* **Changed**
  * [**Issue 142**](https://github.com/tanghaibao/goatools/issues/142)
    * Write GO hierachy to a file now writes a file when using Python3
    * Gaf reader defaults to gaf file version of 2.1 if no version line if found

Release 2019-09-29 0.9.9
-------------------------

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
* **Added**`
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
  * Aspect counts (BP, MF, CC totals) in *TermCounts* object

Release 2019-07-26 0.9.7
-------------------------

### Summary
* **Added**    
  * Support traversing optional relationships when propagating counts for GOEAs

### Details

- Support traversing optional relationships,
  like *part_of* and *regulates*, when doing propagating counts
  using these two new arguments in scripts/find_enrichment.py:
```
        -r, --relationship    Propagate counts up all relationships (default: False)
       --relationships        [RELATIONSHIPS [RELATIONSHIPS ...]]
                               Propagate counts up user-specified relationships
                               (default: None)
```
https://github.com/tanghaibao/goatools/issues/117#issuecomment-515484624    

- Find all ancesters of a GO term using a user-specified list of relationships    
  https://github.com/tanghaibao/goatools/issues/126#issuecomment-511985524    


Release 2019-05-08 0.9.5
-------------------------

### Summary
* **Added**    
  * [**Issue 119: Added support for specifying specific evidence codes**](#added-support-for-specifying-specific-evidence-codes)    
  * [**Issue 127: Added splitting GOEA into BP, MF, and CC**](#added-splitting-goea-into-bp-mf-and-cc)    
* **Fixed**    
  * [**Issue 120 and 121: Lin similarity is positive**](https://github.com/tanghaibao/goatools/issues/120)

### Details
#### Added support for specifying specific evidence codes:
- Specify evidence codes, like EXP (Inferred from Experiment), to exclude or include in a GOEA.
- Specify evidence classes, like Experimental (EXP IDA IPI IMP IGI IEP), which include many evidence codes.
- Get evidence code help:
```
  $ python3 scripts/find_enrichment.py --ev_help
  $ python3 scripts/find_enrichment.py --ev_help_short
```
https://github.com/tanghaibao/goatools/issues/119#issuecomment-488508518    
https://github.com/tanghaibao/goatools/issues/119#issuecomment-489816413    

#### Added splitting GOEA into BP, MF, and CC
Split GOEA into three separate analyses by default:
  * BP (biological process)
  * MF (molecular function)
  * CC cellular component)    

https://github.com/tanghaibao/goatools/issues/127#issuecomment-489776548    


Release 2018-11-27 0.8.12
-------------------------

Release 2018-11-21 0.8.11
-------------------------

Release 2018-09-09 0.8.9
-------------------------

Release 2018-04-28 0.8.4
-------------------------

Release 2018-02-22 0.8.2
-------------------------

Release 2018-02-19 0.7.11
-------------------------

Release 2018-01-07 0.7.11
-------------------------

Release 2017-10-14 0.7.11
-------------------------

Release 2017-10-07 0.7.10
-------------------------

Release 2015-02-05
-------------------------

Release 2014-08-10
-------------------------

Release 2010-04-12
-------------------------

Release 2010-02-23
-------------------------


git log --since=2019-05-08 --before=2019-07-26
