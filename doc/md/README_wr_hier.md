# Print all or part of the GO hierarchy

  * [**Write the entire GO hierarchy for all genes from one species (human)**](#write-the-entire-go-hierarchy-for-all-genes-from-one-species-human)


## Write the entire GO hierarchy for all genes from one species (human)

```
$ scripts/wr_hier.py BP MF CC --gene2go=gene2go --taxid=9606 --dash_len=17 --concise -o human_BP_MF_CC.txt
```
### Arguments

  * **BP MF CC** => Aliases for top GO IDs for biological_process(GO:0008150), molecular_function(GO:0003674) & cellular_process(GO:0005575)

  * **--gene2go=gene2go --taxid=9606** => Tells wr_hier to read the associations in file, "gene2go" for human(9606)

  * **--dash_len=17** => Adds spacing so the GO information printed next to each GO ID does not shift around. The BP branch is 17 GO IDs deep, so 17 is the right number for this print.

  * **--concise** => Because GO IDs can have multiple parents, the same sub-path can be printed multiple times, taking up gobs of space in a report and making the hierarchy report more difficult to read. To print a sub-path just once when it occurs more than one time, use this arg. If you use this arg, your report will be 37k lines. If you do not use this arg, your report will be 351,557 lines due to duplication

### Further Details

GO IDs that are in the association will be marked with a ">" at the beginning of
the line. The numbers on the GO line include the total count of all descendants
(e.g., 265) determined from the GO DAG and the information content score
determined from the associations (e.g., 7.01)

```
  --- GO:0009628   BP 265  7.01 D02 response to abiotic stimulus
  ---- GO:0071214  BP  68  8.30 D04 cellular response to abiotic stimulus
> ----- GO:0071257 BP   0 11.48 D05 cellular response to electrical stimulus
```

The two ancestor GO IDs preceding the marked(>) GO:0071257 are not in the human
association, but are in the hierarchy leading down to GO:0071257 , which is in
the association.

**--concise** causes a line which uses '=' to indicate depth instead of '-':

```
> ======== GO:0071295  BP 13 10.97 L05 D07 cellular response to vitamin
```

The use of '=' instead of '-' to show a level of depth means that this GO AND
its lower level terms were printed earlier and will not be duplicated to save
space. In the case of GO:0071295, the first detailed print looks like this:

```
> -------- GO:0071295  BP 13 10.97 L05 D07 cellular response to vitamin
> --------- GO:0071307 BP  0 13.27 L06 D08 cellular response to vitamin K
> --------- GO:0071231 BP  0 13.27 L05 D08 cellular response to folic acid
> --------- GO:0071301 BP  0 13.96 L05 D08 cellular response to vitamin B1
> --------- GO:0071298 BP  0 13.96 L05 D08 cellular response to L-ascorbic acid
> --------- GO:0071305 BP  0 11.56 L05 D08 cellular response to vitamin D
> --------- GO:0071306 BP  0 13.27 L05 D08 cellular response to vitamin E
```

Copyright (C) 2010-2018, DV Klopfenstein, Haibao Tang et al. All rights reserved.
