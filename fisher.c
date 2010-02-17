/*
 Statistical assessment of an over- or under-representation of a property
 in a list of objects using Fisher's exact test  */

#include <stdio.h>
#include <math.h>

#define min(a,b) ((a)<(b))?(a):(b)
#define max(a,b) ((a)>(b))?(a):(b)

static double __hypergeometric_probability (int i, int n, int C, int G);
static double __lncombination (int n, int p);
static double  __lnfactorial (int n);
static double __lngamma (int z);

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*
	# Compute the hypergeometric probability, or probability that a list of
	# 'n' objects should contain 'i' ones with a particular property when the
	# list has been selected randomly without replacement from a set of 'G'
	# objects in which 'C' exhibit the same property
*/
static double __hypergeometric_probability (int i, int n, int C, int G)
{
    return exp(
      __lncombination(C, i) +
      __lncombination(G - C, n - i) -
      __lncombination(G, n)
     );
}

/*
	# Logarithm of the number of combinations of 'n' objects taken 'p' at a time
*/
static double __lncombination (int n, int p)
{
    return 
      __lnfactorial(n) - 
      __lnfactorial(p) - 
      __lnfactorial(n - p);
}

/* 
	# Logarithm of n! with algorithmic approximation
	# Reference:
	#   Lanczos, C. 'A precision approximation of the gamma function',
	#   J. SIAM Numer. Anal., B, 1, 86-96, 1964."
	#   http://www.matforsk.no/ola/fisher.htm 
*/
static double __lnfactorial (int n)
{
    return (n <= 1)?0: __lngamma(n+1);
}

static double __lngamma (int z)
{
    double x = 0;
    x += 0.1659470187408462e-06 / (z + 7);
    x += 0.9934937113930748e-05 / (z + 6);
    x -= 0.1385710331296526 / (z + 5);
    x += 12.50734324009056 / (z + 4);
    x -= 176.6150291498386 / (z + 3);
    x += 771.3234287757674 / (z + 2);
    x -= 1259.139216722289 / (z + 1);
    x += 676.5203681218835 / (z);
    x += 0.9999999999995183;

    return log(x) - 5.58106146679532777 - z + (z - 0.5) * log(z + 6.5);
}

/*
	# Calculate the Fisher's exact test.
	# Return the left, right and two-tailed p-values
	#
	# Adapted from WordHoard project - http://wordhoard.northwestern.edu
	# (edu.northwestern.at.utils.math.statistics.FishersExactTest)
*/
void pvalue (int k, int n, int C, int G, 
        double *left_tail, double *right_tail, double *two_tailed)
{
    int um = min(n, C), lm = max(0, n+C-G);

    if (um == lm)
    {
        *left_tail = *right_tail = *two_tailed = 1.0;
        return;
    }

	double cutoff = __hypergeometric_probability(k, n, C, G);
	*left_tail = *right_tail = *two_tailed = 0;

    int i;
    for (i=lm; i<um+1; i++)
    {
        double p = __hypergeometric_probability(i, n, C, G);

        if (i <= k) *left_tail += p;
        if (i >= k) *right_tail += p;
        if (p <= cutoff) *two_tailed += p;
    }

    *left_tail = min(*left_tail, 1);
    *right_tail = min(*right_tail, 1);
    *two_tailed = min(*two_tailed, 1);

    return;
}

/*  
	# Enrichment for the attribute in the query
	# (< 1 if the query is impoverished)
*/
double enrichment (int k, int n, int C, int G)
{
    return ((double) k / C) / ((double) n / G);
}


/* print tabular summary of result */
void print_summary (int k, int n, int C, int G)
{
    printf(" Fisher's exact test\n\n");
    printf("  Contingency table :\n\n");
    printf("              Having the | Not having |\n");
    printf("               property  |            | Total\n");
    printf("     Selected   %5d (k)    %5d        %5d (n)\n", k, n - k, n);
    printf(" Not selected   %5d        %5d        %5d\n", C - k, G - C - n + k, G - n);
    printf("        Total   %5d (C)    %5d        %5d (G)\n\n", C, G - C, G);
    printf("         Enrichment : %.6lf", enrichment(k, n, C, G));

    double left, right, two_tailed;
    pvalue(k, n, C, G, &left, &right, &two_tailed);

    printf("       Left p-value : %.6g\n", left);
    printf("      Right p-value : %.6g\n", right);
    printf(" Two-tailed p-value : %.6g\n", two_tailed);

    return;
}

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

int main (int argc, char const* argv[])
{
    print_summary(10, 10, 20, 100);
    
    return 0;
}
