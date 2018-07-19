# Print all or part of the GO hierarchy

  * [**Write the entire GO hierarchy for all genes from one species (human)**]


## Write the entire GO hierarchy for all genes from one species (human)

```
$ scripts/wr_hier.py BP MF CC --gene2go=gene2go --taxid=9606 --dash_len=17 --concise -o human_BP_MF_CC.txt
```

  * **BP MF CC** => Aliases for top GO IDs for biological_process(GO:0008150), molecular_function(GO:0003674) & cellular_process(GO:0005575)

  * **BP MF CC** => Aliases for top GO IDs for biological_process(GO:0008150), molecular_function(GO:0003674) & cellular_process(GO:0005575)

  * **--gene2go=gene2go --taxid=9606** => Tells wr_hier to read the associations in file, "gene2go" for human(9606)

  * **--dash_len=17** => Adds spacing so the GO information printed next to each GO ID does not shift around. The BP branch is 17 GO IDs deep, so 17 is the right number for this print.

  * **--concise** => Because GO IDs can have multiple parents, the same sub-path can be printed multiple times, taking up gobs of space in a report and making the hierarchy report more difficult to read. To print a sub-path just once when it occurs more than one time, use this arg. If you use this arg, your report will be 37k lines. If you do not use this arg, your report will be 351,557 lines due to duplication


Copyright (C) 2010-2018, DV Klopfenstein, Haibao Tang et al. All rights reserved.
