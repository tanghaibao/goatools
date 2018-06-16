# Plot the GO hierarchy with go_plot.py

  * [Plot one GO term and its ancestors]()
  * [Plot one GO term and its ancestors (with relationships)]()

## Plot one GO term and its ancestors
Plot one term and all ancestors using the 'is_a' attribute.    
```
scripts/go_plot.py GO:0003304
```
![heart_jogging_r0](../images/plot_go/GO_0003304_myocardial_epithelial_involution_involved_in_heart_jogging.png)


## Plot one GO term and its ancestors (with relationships)
Plot one term and all ancestors using the 'is_a' attribute and all relationships (--r).

```
scripts/go_plot.py GO:0003304 --r
```
![heart_jogging_r1](../images/plot_go/GO_0003304_myocardial_epithelial_involution_involved_in_heart_jogging_r1.png)

Copyright (C) 2010-2018, DV Klopfenstein, Haibao Tang et al. All rights reserved.
