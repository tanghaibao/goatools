%module fisher
%include "typemaps.i"

void print_summary(int k,int n,int C,int G);
double enrichment(int k,int n,int C,int G);
void pvalue(int k,int n,int C,int G,double *OUTPUT,double *OUTPUT,double *OUTPUT);
