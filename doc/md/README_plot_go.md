# Plot the GO hierarchy with go_plot.py

  * Plot one GO term and its ancestors:
    * [Using 'is_a' relationship only (default)](#plot-one-go-term-and-its-ancestors)
    * [Using all relationships](#plot-one-go-term-and-its-ancestors-with-relationships)
  * Plot two GO terms and their ancestors:
    * [Using 'is_a' relationship only (default)](#plot-two-go-terms-and-their-ancestors)
    * [Using all relationships](#plot-two-go-terms-and-their-ancestors-with-relationships)

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


## Plot two GO terms and their ancestors
Plot two terms and all ancestors using the 'is_a' attribute.        
```
scripts/go_plot.py GO:0003304 GO:0003146 -o heart_jogging.png
```

![heart_jogging_r0](../images/plot_go/heart_jogging.png)


## Plot two GO terms and their ancestors (with relationships)
Plot two terms and all ancestors using the 'is_a' attribute and all relationships (--r).

```
scripts/go_plot.py GO:0003304 GO:0003146 --r -o heart_jogging_r1.png
```
![heart_jogging_r1](../images/plot_go/heart_jogging_r1.png


Copyright (C) 2010-2018, DV Klopfenstein, Haibao Tang et al. All rights reserved.
