# NCBI GeneID and Symbol
NCBI Gene data changes daily, so download gene lists frequently.    

## Example of printing gene symbols (instead of geneids) from NCBI gene in a GOEA results table.

1. [Download **NCBI gene** information into **gene_result.txt**](#1-download-ncbi-gene-information-into-gene_resulttxt)
2. [Convert **gene_result.txt** file into Python](#2-convert-gene_resulttxt-file-into-python)
3. [Create a geneid2symbol.txt ASCII file](#3-create-a-geneid2symboltxt-ascii-file)
4. [Run GOEA using **--id2name** flag](#4-run-goea-using---id2name-flag)

## 1) Download NCBI gene information into gene_result.txt

  * 1a) [Search for a gene set](#1a-search-for-a-gene-set)    
  * 1b) [Download gene set into gene_result.txt](#1b-download-gene_resulttxt)    

### 1a) Search for a gene set

 * Go to NCBI Gene in your browser https://www.ncbi.nlm.nih.gov/gene/
 * Type a search pattern and hit the "Search" button.    
   For example, to download all human protein-coding genes:
```
Text in 'Search':
genetype protein coding[Properties] AND "9606"[Taxonomy ID] AND alive[property] 
```
![NCBI Gene Search](/doc/images/NCBI_gene_search.png)

    | NCBI Search text                    | Description
    |-------------------------------------|----------------------
    | genetype protein coding[Properties] | Protein-coding genes
    | "9606"[Taxonomy ID]                 | human
    | alive[property]                     | NOT obsolete

### 1b) Download gene_result.txt
  * Click pull-down menu: "**Send to:**"    
  * Click "**File**" radial button    
  * Click "**Create File**" button    
![NCBI Gene Download](/doc/images/NCBI_gene_download.png)


## 2) Convert gene_result.txt file into Python 
```
$ scripts/ncbi_gene_results_to_python.py -i gene_result.txt -o gene_result.py
      20380 lines READ:  gene_result.txt
      20361 geneids WROTE: gene_result.py
```

## 3) Create a geneid2symbol.txt ASCII file
COMING SOON...
## 4) Run GOEA using --id2name flag
COMING SOON...


Copyright (C) 2016-2018, DV Klopfenstein, Haibao Tang et al. All rights reserved.
