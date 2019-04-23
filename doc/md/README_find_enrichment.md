# Gene Ontolog Enrichment Analyses (GOEA) from the command-line

## Run GOEA and print GO terms

  1. with pvalues < 0.05 corrected using Benjamini-Hochberg multiple test correction using annotation in a:    
     * Original id-to-GOs text format   
     * GAF file    
     * GPAD file    
     * NCBI's gene2go file
  2. [with **uncorrected pvalues < 0.05** to the screen](#2-print-go-terms-with-uncorrected-pvalues--005-to-the-screen-default):
    (default)
  3. [with **uncorrected pvalues < 0.05** to the screen, **grouped**](#3-print-go-terms-with-uncorrected-pvalues--005-to-the-screen-grouped):    
     **--sections=...**
  4. [with **uncorrected pvalues < 0.05** to files](#4-print-go-terms-with-uncorrected-pvalues--005-to-an-xlsx-file):    
    **--outfile=goea.xlsx**
  5. [with **Benjamini-Hochberg pvalues < 0.05** (plus bonferroni, sidak, holm)](#5-print-go-terms-with-benjamini-hochberg-pvalues--005) 
    **--pval_field=fdr_bh**    
  6. [with **Benjamini-Hochberg pvalues < 0.05** to an **xlsx** file](#6-print-go-terms-with-benjamini-hochberg-only-pvalues--005-to-an-xlsx-file)     
    **--method=fdr_bh --outfile=goea_fdr_bh_flat.xlsx**    
  7. [with **Benjamini-Hochberg pvalues < 0.05** to an **xlsx** file, **grouped**](#7-print-go-terms-with-benjamini-hochberg-only-pvalues--005-to-an-xlsx-file-grouped)     
    **--method=fdr_bh --outfile=goea_fdr_bh_grpd.xlsx --sections...**    
  8. [regardless of pvalue **(ALL GO terms)**](#8-print-all-go-terms-regardless-of-pvalue)
    **--pval=-1**    

## Example Details

CMD: python scripts/find_enrichment.py data/study data/population data/association --outfile=goea.xlsx,goea.tsv --pval_field=fdr_bh     

### 1) Enrichment Analysis using various annotation formats    

Print results where Benjamini-Hochberg values are less than 0.05    
```
--pval=0.05             # print pvalue's < 0.05
--method=fdr_bh         # Benjamini-Hochberg multiple test correction on uncorrected p-values
--pval_field=fdr_bh     # print fdr-bh values < 0.05 (rather than uncorrected pvalues)
--outfile=FILENAME.xlsx # Write to an Excel spreadsheet
```

### 1A) Original id-to-GOs text format    

### 1B) GAF file

**COMMAND:**    
```
python3 scripts/find_enrichment.py ids_stu_gaf.txt ids_pop_gaf.txt goa_human.gaf --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --outfile=results_gaf.xlsx
```

**RESULTS:**    
```
go-basic.obo: fmt(1.2) rel(2019-04-17) 47,398 GO Terms
HMS:0:00:13.741054 424,966 annotations READ: goa_human.gaf
Study: 100 vs. Population 19427

fisher module not installed.  Falling back on scipy.stats.fisher_exact
Propagating term counts to parents ..
1 GO IDs in assc. are not found in the GO-DAG: GO:0034437
100% 19,366 of 19,427 population items found in association
 94%     94 of    100 study items found in association
100%    100 of    100 study items found in population(19427)
Calculating 21,575 uncorrected p-values using fisher_scipy_stats
  21,575 GO terms are associated with 19,366 of 19,427 population items
   2,055 GO terms are associated with     94 of    100 study items
       7 GO terms found significant (< 0.05=alpha) (  5 enriched +   2 purified): statsmodels fdr_bh
      49 study items associated with significant GO IDs (enriched)
      19 study items associated with significant GO IDs (purified)
    481 of 21,575 results have uncorrected P-values <= 0.05=pval

      7 items WROTE: results_gaf.xlsx
```


### 1B) GPAD file

### 1D) NCBI's gene2go file

### 2) Print GO terms with uncorrected pvalues < 0.05 to the screen (default)
The default is to print all GO terms with uncorrected P-values < 0.05 (default)

```
python scripts/find_enrichment.py data/study data/population data/association
```


|GO|NS|ep|name|study|pop|p_uncorrected|depth|# study|p_bonferroni|p_sidak|p_holm|p_fdr_bh|study_items
|--|--|--|----|-----|---|-------------|-----|-----------|------------|-------|------|--------|-----------
|GO:0036211|BP|e|protein modification process|33/276|1725/33239|8.070e-06|5|33|0.0491|0.04789|0.04910|0.00982|AT1G13580...
|GO:0006464|BP|e|cellular protein modification process|33/276|1725/33239|8.070e-06|6|33|0.0491|0.04789|0.04910|0.00982|AT1G13580...

### 3) Print GO terms with uncorrected pvalues < 0.05 to the screen, grouped 
optional attribute: **--sections=goatools.test_data.sections.data2018_07_find_enrichment**

Acceptable values for **--sections** are:

| File type     | Example value
|---------------|------------------------------------------------------------
| Python module | goatools.test_data.sections.data2018_07_find_enrichment
| text file     | sections.txt
| Python file   | goatools/test_data/sections/data2018_07_find_enrichment.py


```
python scripts/find_enrichment.py data/study data/population data/association --sections=goatools.test_data.sections.data2018_07_find_enrichment

protein
GO:0036211 BP e 9.83e-03  1320  5.14 D05  33/276  1725/33239 protein modification process ...
GO:0006464 BP e 9.83e-03  1318  5.14 D06  33/276  1725/33239 cellular protein modification process ...
...

phosph
GO:0016301 MF e 4.64e-02   355  4.87 D04  25/276  1310/33239 kinase activity ...
...

```

### 4) Print GO terms with uncorrected pvalues < 0.05 to an xlsx file
optional attribute: **--outfile=goea_uncorr.xlsx**    

```
$ python scripts/find_enrichment.py data/study data/population data/association --outfile=goea_uncorr.xlsx
    ...
    253 items WROTE: goea.xlsx
```

### 5) Print GO terms with Benjamini-Hochberg pvalues < 0.05
optional attribute: **--pval_field=fdr_bh**     
optional attribute: **----outfile=goea_fdr_bh.xlsx,goea_fdr_bh.tsv**    

```
$ python scripts/find_enrichment.py data/study data/population data/association --outfile=goea_fdr_bh.xlsx,goea_fdr_bh.tsv --pval_field=fdr_bh
     17 items WROTE: goea.xlsx
     17 items WROTE: goea.tsv
```

### 6) Print GO terms with Benjamini-Hochberg (only) pvalues < 0.05 to an xlsx file
optional attribute: **--method=fdr_bh**    
optional attribute: **--outfile=goea_fdr_bh_flat.xlsx**    

```
$ python scripts/find_enrichment.py data/study data/population data/association --outfile=goea_fdr_bh_flat.xlsx --method=fdr_bh
     17 items WROTE: goea_fdr_bh_flat.xlsx
```

### 7) Print GO terms with Benjamini-Hochberg (only) pvalues < 0.05 to an xlsx file, grouped
optional attribute: **--method=fdr_bh**    
optional attribute: **--outfile=goea_fdr_bh_grpd.xlsx**    
optional attribute: **--sections=goatools.test_data.sections.data2018_07_find_enrichment**    

```
$ python scripts/find_enrichment.py data/study data/population data/association --outfile=goea_fdr_bh_grpd.xlsx --method=fdr_bh --sections=goatools.test_data.sections.data2018_07_find_enrichment
     17 items WROTE: goea_fdr_bh_grpd.xlsx
```

### 8) Print all GO terms, regardless of pvalue
optional attribute: **--pval=-1**

```
$ python scripts/find_enrichment.py data/study data/population data/association --outfile=goea_all.xlsx,goea_all.tsv --pval=-1
   6088 items WROTE: goea.xlsx
   6088 items WROTE: goea.tsv
```

Copyright (C) 2010-2018, DV Klopfenstein, Haibao Tang, et al. All rights reserved.
