# Gene Ontolog Enrichment Analyses (GOEA) from the command-line

## GOEA command-line examples

  1. [with **pvalues < 0.05 corrected using Benjamini-Hochberg multiple test
		 correction**](#1-enrichment-analysis-using-various-annotation-formats) using annotation in a:    
     * Original id-to-GOs text format   
     * GAF file    
     * GPAD file    
     * NCBI's gene2go file
  2. **Exclude or include annotations by evidence code**    
     Use with GAF, GPAD, or NCBI's gene2go files:
     * [**Exclude annotations inferred from Electronic Annotation (IEA)**](#2a-exclude-annotations-inferred-from-electronic-annotation-iea)
     * [**Include only annotations inferred from experimental evidence**](#2b-include-only-annotations-inferred-from-experimental-evidence)
     * [Get **list of evidence codes**](#2c-Get-list-of-evidence-codes)
  3. [**Limit the GOEA to _biological process_, _molecular function_, and/or _cellular compartment_**](https://github.com/tanghaibao/goatools/blob/master/doc/md/README_find_enrichment.md#3-limit-the-goea-to-biological-process-molecular-function-andor-cellular-compartment)
  4. [with **uncorrected pvalues < 0.05** to the screen](#4-print-go-terms-with-uncorrected-pvalues--005-to-the-screen-default):
    (default)
  5. [with **uncorrected pvalues < 0.05** to the screen, **grouped**](#5-print-go-terms-with-uncorrected-pvalues--005-to-the-screen-grouped):    
     **--sections=...**
  6. [with **uncorrected pvalues < 0.05** to files](#6-print-go-terms-with-uncorrected-pvalues--005-to-an-xlsx-file):    
    **--outfile=goea.xlsx**
  7. [with **Benjamini-Hochberg pvalues < 0.05** (plus bonferroni, sidak, holm)](#7-print-go-terms-with-benjamini-hochberg-pvalues--005) 
    **--pval_field=fdr_bh**    
  8. [with **Benjamini-Hochberg pvalues < 0.05** to an **xlsx** file](#8-print-go-terms-with-benjamini-hochberg-only-pvalues--005-to-an-xlsx-file)     
    **--method=fdr_bh --outfile=goea_fdr_bh_flat.xlsx**    
  9. [with **Benjamini-Hochberg pvalues < 0.05** to an **xlsx** file, **grouped**](#9-print-go-terms-with-benjamini-hochberg-only-pvalues--005-to-an-xlsx-file-grouped)     
    **--method=fdr_bh --outfile=goea_fdr_bh_grpd.xlsx --sections...**    
  10. [regardless of pvalue **(ALL GO terms)**](#10-print-all-go-terms-regardless-of-pvalue)
    **--pval=-1**    

## Example Details

CMD: python scripts/find_enrichment.py data/study data/population data/association --outfile=goea.xlsx,goea.tsv --pval_field=fdr_bh     

### 1) Enrichment Analysis using various annotation formats    

Arguments to print results where Benjamini-Hochberg values are less than 0.05:    
```
--pval=0.05             # print pvalue's < 0.05
--method=fdr_bh         # Use Benjamini-Hochberg multiple test correction on uncorrected p-values
--pval_field=fdr_bh     # print fdr_bh values < 0.05 (rather than uncorrected pvalues)
--outfile=FILENAME.xlsx # Write to an Excel spreadsheet
```

### 1A) Original id-to-GOs text format    
```
python3 scripts/find_enrichment.py data/study.txt data/population.txt data/association.txt --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --outfile=results_id2gos.xlsx
```

### 1B) GAF file
```
python3 scripts/find_enrichment.py ids_stu.txt ids_pop.txt goa_human.gaf --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --outfile=results_gaf.xlsx
```

### 1B) GPAD file
```
python3 scripts/find_enrichment.py ids_stu.txt ids_pop.txt goa_human.gpad --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --outfile=results_gpad.xlsx
```

### 1D) NCBI's gene2go file
#### Human is the default:
```
python3 scripts/find_enrichment.py ids_stu_gene2go_9606.txt ids_pop_gene2go_9606.txt gene2go --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --outfile=results_gene2go_9606.xlsx
```

#### Specify mouse using the **--taxid** argment:
```
python3 scripts/find_enrichment.py ids_stu_gene2go_10090.txt ids_pop_gene2go_10090.txt gene2go --taxid=10090 --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --outfile=results_gene2go_10090.xlsx
```

### 2) Exclude or include annotations by evidence code

#### 2a) Exclude annotations inferred from Electronic Annotation (IEA)
```
EXCLUDE all IEA annotations:

  --ev_exc=IEA  # Evidence codes with IEA are excluded
```

```
python3 scripts/find_enrichment.py ids_stu_gpad.txt ids_pop_gpad.txt goa_human.gpad --ev_exc=IEA --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --outfile=results_gpad.xlsx
```

#### 2b) Include only annotations inferred from experimental evidence
```
INCLUDE evidence codes by group of by code:

  1) --ev_inc=Experimental             # By group name
  1) --ev_inc=Experimental,Similarity  # By group names
  2) --ev_inc=EXP,IDA,IPI,IMP,IGI,IEP  # List all the Experimental codes
```
**Specify the experimental evidence codes by their group name:**
```
python3 scripts/find_enrichment.py ids_stu_gpad.txt ids_pop_gpad.txt goa_human.gpad --ev_inc=Experimental --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --outfile=results_gpad.xlsx
```

**Specify the experimental evidence codes explicitly:**    
```
python3 scripts/find_enrichment.py ids_stu_gpad.txt ids_pop_gpad.txt goa_human.gpad --ev_inc=EXP,IDA,IPI,IMP,IGI,IEP --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --outfile=results_gpad.xlsx
```

#### 2c) Get list of evidence codes
The argument, **--ev_help**, will cause the Evidence codes to be printed:
```
python3 scripts/find_enrichment.py --ev_help
```
##### Evidence Code Help
You can use either group names, like Experimental, or codes on the commandline:    
```
EVIDENCE CODES:
    Experimental:
        EXP Inferred from Experiment
        IDA Inferred from Direct Assay
        IPI Inferred from Physical Interaction
        IMP Inferred from Mutant Phenotype
        IGI Inferred from Genetic Interaction
        IEP Inferred from Expression Pattern
    Similarity:
        ISS Inferred from Sequence or structural Similarity
        ISO Inferred from Sequence Orthology
        ISA Inferred from Sequence Alignment
        ISM Inferred from Sequence Model used in manual assertion
        IGC Inferred from Genomic Context
        IBA Inferred from Biological aspect of Ancestor
        IBD Inferred from Biological aspect of Descendant
        IKR Inferred from phylogenetic determination of loss of key residues (manual assertion)
        IRD Inferred from Rapid Divergence from ancestral sequence (manual assertion)
        IMR Phylogenetic determination of loss of key residues in manual assertion
    Combinatorial:
        RCA Inferred from Reviewed Computational Analysis
    High_Throughput:
        HTP Inferred from High Throughput Experimental
        HDA Inferred from High Throughput Direct Assay
        HMP Inferred from High Throughput Mutant Phenotype
        HGI Inferred from High Throughput Genetic Interaction
        HEP Inferred from High Throughput Expression Pattern
    Author:
        TAS Traceable Author Statement used in manual assertion
        NAS Non-traceable Author Statement used in manual assertion
    Curatorial:
         IC Inferred by Curator
    No biological data:
         ND No biological Data available
    Automatic:
        IEA Inferred from Electronic Annotation
```

### 3) Limit the GOEA to _biological process_, _molecular function_, and/or _cellular compartment_
Use the --ns option to limit the GOEA to any of:

| NS | Namespace |
|----|--------------|
| BP | Biological Process |
| MF | Molecular Function |
| CC | Cellular Component |

Use the two-letter abbreviation as the value in the --ns option:

```
Namespace examples:
--ns=BP      # Only run one GOEA on the Biological Process branch
--ns=BP,MF   # Run two GOEAs: One on Biological Process and one on MOlecular Function
--ns=CC      # Only run one GOEA on the Molecular Function branch
```

#### Without --ns, all three branches will be analyzed:
```
python3 scripts/find_enrichment.py ids_stu_gene2go_9606.txt ids_pop_gene2go_9606.txt gene2go --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --outfile=results_gene2go_9606.xlsx
```

#### With **--ns=MF**, only the _**molecular function**_ branch will be analyzed:
```
python3 scripts/find_enrichment.py ids_stu_gene2go_9606.txt ids_pop_gene2go_9606.txt gene2go --ns=MF --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --outfile=results_gene2go_9606.xlsx
```

#### With **--ns=BP,MF**, run GOEAs on both the _**biological process**_ and the _**molecular function**_ branches:
```
python3 scripts/find_enrichment.py ids_stu_gene2go_9606.txt ids_pop_gene2go_9606.txt gene2go --ns=BP,MF --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --outfile=results_gene2go_9606.xlsx
```

#### With **--ns=MF** and **--ev_inc=IPI**: 
Run a GOEA on the _**molecular function**_ branch for just the evidence code, **IPI**.     
**IPI** includes evidence codes which are inferred from Physical Interaction    
```
python3 scripts/find_enrichment.py ids_stu_gene2go_9606.txt ids_pop_gene2go_9606.txt -ns=MF --ev_inc=IPI gene2go --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --outfile=results_gene2go_9606.xlsx
```



### 4) Print GO terms with uncorrected pvalues < 0.05 to the screen (default)
The default is to print all GO terms with uncorrected P-values < 0.05 (default)

```
python scripts/find_enrichment.py data/study.txt data/population.txt data/association.txt
```


|GO|NS|ep|name|study|pop|p_uncorrected|depth|# study|p_bonferroni|p_sidak|p_holm|p_fdr_bh|study_items
|--|--|--|----|-----|---|-------------|-----|-----------|------------|-------|------|--------|-----------
|GO:0036211|BP|e|protein modification process|33/276|1725/33239|8.070e-06|5|33|0.0491|0.04789|0.04910|0.00982|AT1G13580...
|GO:0006464|BP|e|cellular protein modification process|33/276|1725/33239|8.070e-06|6|33|0.0491|0.04789|0.04910|0.00982|AT1G13580...

### 5) Print GO terms with uncorrected pvalues < 0.05 to the screen, grouped 
optional attribute: **--sections=goatools.test_data.sections.data2018_07_find_enrichment**

Acceptable values for **--sections** are:

| File type     | Example value
|---------------|------------------------------------------------------------
| Python module | goatools.test_data.sections.data2018_07_find_enrichment
| text file     | sections.txt
| Python file   | goatools/test_data/sections/data2018_07_find_enrichment.py


```
python scripts/find_enrichment.py data/study.txt data/population.txt data/association.txt --sections=goatools.test_data.sections.data2018_07_find_enrichment

protein
GO:0036211 BP e 9.83e-03  1320  5.14 D05  33/276  1725/33239 protein modification process ...
GO:0006464 BP e 9.83e-03  1318  5.14 D06  33/276  1725/33239 cellular protein modification process ...
...

phosph
GO:0016301 MF e 4.64e-02   355  4.87 D04  25/276  1310/33239 kinase activity ...
...

```

### 6) Print GO terms with uncorrected pvalues < 0.05 to an xlsx file
optional attribute: **--outfile=goea_uncorr.xlsx**    

```
$ python scripts/find_enrichment.py data/study.txt data/population.txt data/association.txt --outfile=goea_uncorr.xlsx
    ...
    253 items WROTE: goea.xlsx
```

### 7) Print GO terms with Benjamini-Hochberg pvalues < 0.05
optional attribute: **--pval_field=fdr_bh**     
optional attribute: **----outfile=goea_fdr_bh.xlsx,goea_fdr_bh.tsv**    

```
$ python scripts/find_enrichment.py data/study.txt data/population.txt data/association.txt --outfile=goea_fdr_bh.xlsx,goea_fdr_bh.tsv --pval_field=fdr_bh
     17 items WROTE: goea.xlsx
     17 items WROTE: goea.tsv
```

### 8) Print GO terms with Benjamini-Hochberg (only) pvalues < 0.05 to an xlsx file
optional attribute: **--method=fdr_bh**    
optional attribute: **--outfile=goea_fdr_bh_flat.xlsx**    

```
$ python scripts/find_enrichment.py data/study.txt data/population.txt data/association.txt --outfile=goea_fdr_bh_flat.xlsx --method=fdr_bh
     17 items WROTE: goea_fdr_bh_flat.xlsx
```

### 9) Print GO terms with Benjamini-Hochberg (only) pvalues < 0.05 to an xlsx file, grouped
optional attribute: **--method=fdr_bh**    
optional attribute: **--outfile=goea_fdr_bh_grpd.xlsx**    
optional attribute: **--sections=goatools.test_data.sections.data2018_07_find_enrichment**    

```
$ python scripts/find_enrichment.py data/study.txt data/population.txt data/association.txt --outfile=goea_fdr_bh_grpd.xlsx --method=fdr_bh --sections=goatools.test_data.sections.data2018_07_find_enrichment
     17 items WROTE: goea_fdr_bh_grpd.xlsx
```

### 10) Print all GO terms, regardless of pvalue
optional attribute: **--pval=-1**

```
$ python scripts/find_enrichment.py data/study.txt data/population.txt data/association.txt --outfile=goea_all.xlsx,goea_all.tsv --pval=-1
   6088 items WROTE: goea.xlsx
   6088 items WROTE: goea.tsv
```

Copyright (C) 2010-present, DV Klopfenstein, Haibao Tang, et al. All rights reserved.
