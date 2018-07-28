# NCBI GeneID and Symbol
Example of printing gene symbols (instead of geneids) from NCBI gene in a GOEA results table.

1. Download NCBI gene information into gene_result.txt
2. Convert gene_result.txt file into Python 
3. Create a geneid2symbol.txt ASCII file
4. Run GOEA using --id2name flag

## Download NCBI gene information into gene_result.txt
 * Go to NCBI Gene in your browser https://www.ncbi.nlm.nih.gov/gene/
 * Type a search pattern and hit the "Search" button.    
   For example, to download all human protein-coding genes:
```
genetype protein coding[Properties] AND "9606"[Taxonomy ID] AND alive[property] 
```

    | NCBI Search text                    | Description
    |-------------------------------------|----------------------
    | genetype protein coding[Properties] | Protein-coding genes
    | "9606"[Taxonomy ID]                 | human
    | alive[property]                     | NOT obsolete

![NCBI Gene Search](doc/images/NCBI_gene_search.png)

## Convert gene_result.txt file into Python 
## Create a geneid2symbol.txt ASCII file
## Run GOEA using --id2name flag


Copyright (C) 2016-2018, DV Klopfenstein, Haibao Tang et al. All rights reserved.
