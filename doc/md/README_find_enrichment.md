# Gene Ontolog Enrichment Analyses (GOEA) from the command-line

## Run GOEA and print GO terms

  1. [with **uncorrected pvalues < 0.05** to the screen](#1-print-go-terms-with-uncorrected-pvalues--005-to-the-screen-default):
    (default)
  2. [with **uncorrected pvalues < 0.05** to the screen, **grouped**](#2-print-go-terms-with-uncorrected-pvalues--005-to-the-screen-grouped):    
     **--sections=...**
  3. [with **uncorrected pvalues < 0.05** to files](#3-print-go-terms-with-uncorrected-pvalues--005-to-an-xlsx-file):    
    **--outfile=goea.xlsx**
  4. [with **Benjamini-Hochberg pvalues < 0.05** (plus bonferroni, sidak, holm)](#4-print-go-terms-with-benjamini-hochberg-pvalues--005) 
    **--pval_field=fdr_bh**    
  5. [with **Benjamini-Hochberg pvalues < 0.05** to an **xlsx** file](#5-print-go-terms-with-benjamini-hochberg-only-pvalues--005-to-an-xlsx-file)     
    **--method=fdr_bh --outfile=goea_fdr_bh_flat.xlsx**    
  6. [with **Benjamini-Hochberg pvalues < 0.05** to an **xlsx** file, **grouped**](#6-print-go-terms-with-benjamini-hochberg-only-pvalues--005-to-an-xlsx-file-grouped)     
    **--method=fdr_bh --outfile=goea_fdr_bh_grpd.xlsx --sections...**    
  7. [regardless of pvalue **(ALL GO terms)**](#7-print-all-go-terms-regardless-of-pvalue)
    **--pval=-1**    

## Example Details

CMD: python scripts/find_enrichment.py data/study data/population data/association --outfile=goea.xlsx,goea.tsv --pval_field=fdr_bh     

### 1) Print GO terms with uncorrected pvalues < 0.05 to the screen (default)
The default is to print all GO terms with uncorrected P-values < 0.05 (default)

```
python scripts/find_enrichment.py data/study data/population data/association
```


|GO|NS|ep|name|study|pop|p_uncorrected|depth|# study|p_bonferroni|p_sidak|p_holm|p_fdr_bh|study_items
|--|--|--|----|-----|---|-------------|-----|-----------|------------|-------|------|--------|-----------
|GO:0036211|BP|e|protein modification process|33/276|1725/33239|8.07048771243e-06|5|33|0.0491331291933|0.0478945023525|0.0491089177301|0.00982662583865|AT1G13580...
|GO:0006464|BP|e|cellular protein modification process|33/276|1725/33239|8.07048771243e-06|6|33|0.0491331291933|0.0478945023525|0.0491089177301|0.00982662583865|AT1G13580...

### 2) Print GO terms with uncorrected pvalues < 0.05 to the screen, grouped 
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

### 3) Print GO terms with uncorrected pvalues < 0.05 to an xlsx file
optional attribute: **--outfile=goea_uncorr.xlsx**    

```
$ python scripts/find_enrichment.py data/study data/population data/association --outfile=goea_uncorr.xlsx
    ...
    253 items WROTE: goea.xlsx
```

### 4) Print GO terms with Benjamini-Hochberg pvalues < 0.05
optional attribute: **--pval_field=fdr_bh**     
optional attribute: **----outfile=goea_fdr_bh.xlsx,goea_fdr_bh.tsv**    

```
$ python scripts/find_enrichment.py data/study data/population data/association --outfile=goea_fdr_bh.xlsx,goea_fdr_bh.tsv --pval_field=fdr_bh
     17 items WROTE: goea.xlsx
     17 items WROTE: goea.tsv
```

### 5) Print GO terms with Benjamini-Hochberg (only) pvalues < 0.05 to an xlsx file
optional attribute: **--method=fdr_bh**    
optional attribute: **--outfile=goea_fdr_bh_flat.xlsx**    

```
$ python scripts/find_enrichment.py data/study data/population data/association --outfile=goea_fdr_bh_flat.xlsx --method=fdr_bh
     17 items WROTE: goea_fdr_bh_flat.xlsx
```

### 6) Print GO terms with Benjamini-Hochberg (only) pvalues < 0.05 to an xlsx file, grouped
optional attribute: **--method=fdr_bh**    
optional attribute: **--outfile=goea_fdr_bh_grpd.xlsx**    
optional attribute: **--sections=goatools.test_data.sections.data2018_07_find_enrichment**    

```
$ python scripts/find_enrichment.py data/study data/population data/association --outfile=goea_fdr_bh_grpd.xlsx --method=fdr_bh --sections=goatools.test_data.sections.data2018_07_find_enrichment
     17 items WROTE: goea_fdr_bh_grpd.xlsx
```

### 7) Print all GO terms, regardless of pvalue
optional attribute: **--pval=-1**

```
$ python scripts/find_enrichment.py data/study data/population data/association --outfile=goea_all.xlsx,goea_all.tsv --pval=-1
   6088 items WROTE: goea.xlsx
   6088 items WROTE: goea.tsv
```

Copyright (C) 2010-2018, DV Klopfenstein, Haibao Tang, et al. All rights reserved.
