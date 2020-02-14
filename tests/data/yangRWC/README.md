# Yang random walk contribution (RWC) tests
[**Improving GO semantic similarity measures by exploring the ontology beneath the terms and modelling uncertainty**](https://pubmed.ncbi.nlm.nih.gov/22522134)    
Yang, Haixuan et al. Bioinformatics (2012)    

Yang's innovative addition to similarity calculation between two GO terms is to
take into account the descendants of the two GO terms,
leading to greatly imporved similarity measures.


## GO DAGs used to test the Python port of Yang RWC
The [**high-quality Java implementation**](https://github.com/pwac092/gossto) of Yang's RWC semantic similarity addition is described in:    
[**GOssTo: a stand-alone application and a web tool for calculating semantic similarities on the Gene Ontology**](https://pubmed.ncbi.nlm.nih.gov/24659104)   
by first authors,
[**Alfonso E Romero**](https://pubmed.ncbi.nlm.nih.gov/?term=Romero+AE%5BAU%5D) and
[**Horacio Caniza**](https://pubmed.ncbi.nlm.nih.gov/?term=Caniza+H%5BAU%5D).

### Original Yang's Fig 2
![Yang Fig2](bioinf_yang_fig2.png)

#### Yang's Fig 2a GO DAG in GOATOOLS
A total of 100 genes are associated with GO terms.

Note that we swapped node F and G compared to Yang Fig 2.    
![fig2a](yang_fig2a.png)

#### Fig2a modified: Reduce number of annotations by 10x
We reduced 100 genes to 10 genes (a, b, c, .... h, i, j) for easier testing.

![fig2a_small](yang_fig2a_small.png)

#### Fig2a modified: Gene, a, is explicitly annotated to both B and D
![fig2a_small_b](yang_fig2a_small_b.png)

#### Fig2a modified: F is not annotated
![fig2a_nonleaf0](yang_fig2a_nonleaf0.png)

#### Simplify: All GO terms have one parent. F is not annotated.
![faa2a_nonleaf0](yang_faa2a_nonleaf0.png)


### Original Yang's Fig 1
![Yang Fig 1](bioinf_yang_fig1.png)

#### Yang Fig 1a in GOATOOLS
![fig1a](yang_fig1a.png)
![fig1a_small](yang_fig1a_small.png)

#### Yang Fig 1b in GOATOOLS
![fig1b](yang_fig1b.png)
