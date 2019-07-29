# CHANGELOG

* [Release 2019-07-26 0.9.7](#release-2019-07-26-097)
* [Release 2019-05-08 0.9.5](#release-2019-05-08-095)

Latest changes
--------------

Release 2019-07-26 0.9.7
-------------------------

### Added:

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

- https://github.com/tanghaibao/goatools/issues/126#issuecomment-511985524    
  Find all ancesters of a GO term using a user-specified list of relationships


Release 2019-05-08 0.9.5
-------------------------

### Changes
* **Added**    
  * [**Issue 119: Added support for specifying specific evidence codes**](#added-support-for-specifying-specific-evidence-codes)    
  * [**Issue 127: Added splitting GOEA into BP, MF, and CC**](#added-splitting-goea-into-bp-mf-and-cc)    
* **Fixed**    
  * [**Issue 120 and 121: Lin similarity is positive**](https://github.com/tanghaibao/goatools/issues/120)

### Added support for specifying specific evidence codes:
https://github.com/tanghaibao/goatools/issues/119#issuecomment-488508518    
https://github.com/tanghaibao/goatools/issues/119#issuecomment-489816413    

- Specify evidence codes, like EXP (Inferred from Experiment), to exclude or include in a GOEA.
- Specify evidence code classes, like Experimental (EXP IDA IPI IMP IGI IEP).
- Get evidence code help:
```
  $ python3 scripts/find_enrichment.py --ev_help
  $ python3 scripts/find_enrichment.py --ev_help_short
```

### Added splitting GOEA into BP, MF, and CC
https://github.com/tanghaibao/goatools/issues/127#issuecomment-489776548    

Split GOEA into three separate analyses by default:
  * BP (biological process)
  * MF (molecular function)
  * CC cellular component)    



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
