
# Statistical assessment of an over- or under-representation of a property
# in a list of objects using Fisher's exact test

from math import log, exp

class FisherExactTest:

	# Calculate the Fisher's exact test.
	# Return the left, right and two-tailed p-values
	#
	# Adapted from WordHoard project - http://wordhoard.northwestern.edu
	# (edu.northwestern.at.utils.math.statistics.FishersExactTest)
	def pvalue (self, k, n, C, G):
		um, lm = min(n, C), max(0, n + C - G)

		if (um == lm):
			return 1.0, 1.0, 1.0

		cutoff = self.__hypergeometric_probability(k, n, C, G)
		left_tail, right_tail, two_tailed = 0, 0, 0

		for i in range(lm, um + 1):
			p = self.__hypergeometric_probability(i, n, C, G)

			if (i <= k):
				left_tail += p

			if (i >= k):
				right_tail += p

			if (p <= cutoff):
				two_tailed += p

		left_tail = min(left_tail, 1)
		right_tail = min(right_tail, 1)
		two_tailed = min(two_tailed, 1)

		return left_tail, right_tail, two_tailed

	# Enrichment for the attribute in the query
	# (< 1 if the query is impoverished)
	def enrichment (self, k, n, C, G):
		return (float(k) / C) / (float(n) / G)

	def evaluate (self, k, n, C, G):
		return self.enrichment(k, n, C, G), self.pvalue(k, n, C, G)

	def print_report (self, k, n, C, G):
		print " Fisher's exact test\n"
		print "  Contingency table :\n"
		print "              Having the | Not having |"
		print "               property  |            | Total"
		print "     Selected   %5s (k)    %5s        %5s (n)" % (k, n - k, n)
		print " Not selected   %5s        %5s        %5s" % (C - k, G - C - n + k, G - n)
		print "        Total   %5s (C)    %5s        %5s (G)" % (C, G - C, G)
		print
		print "         Enrichment : %.6f" % self.enrichment(k, n, C, G)

		left, right, two_tailed = self.pvalue(k, n, C, G)

		print "       Left p-value : %.6g" % left
		print "      Right p-value : %.6g" % right
		print " Two-tailed p-value : %.6g" % two_tailed

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	# Compute the hypergeometric probability, or probability that a list of
	# 'n' objects should contain 'i' ones with a particular property when the
	# list has been selected randomly without replacement from a set of 'G'
	# objects in which 'C' exhibit the same property
	def __hypergeometric_probability (self, i, n, C, G):
		return exp(
		  self.__lncombination(C, i) +
		  self.__lncombination(G - C, n - i) -
		  self.__lncombination(G, n)
		 )

	# Logarithm of the number of combinations of 'n' objects taken 'p' at a time
	def __lncombination (self, n, p):
		return \
		  self.__lnfactorial(n) - \
		  self.__lnfactorial(p) - \
		  self.__lnfactorial(n - p)

	# Logarithm of n! with algorithmic approximation
	# Reference:
	#   Lanczos, C. 'A precision approximation of the gamma function',
	#   J. SIAM Numer. Anal., B, 1, 86-96, 1964."
	#   http://www.matforsk.no/ola/fisher.htm 
	def __lnfactorial (self, n):
		if (n <= 1):
			return 0
		else:
			return self.__lngamma(n + 1)

	def __lngamma (self, z):
		x = 0
		x += 0.1659470187408462e-06 / (z + 7)
		x += 0.9934937113930748e-05 / (z + 6)
		x -= 0.1385710331296526 / (z + 5)
		x += 12.50734324009056 / (z + 4)
		x -= 176.6150291498386 / (z + 3)
		x += 771.3234287757674 / (z + 2)
		x -= 1259.139216722289 / (z + 1)
		x += 676.5203681218835 / (z)
		x += 0.9999999999995183

		return log(x) - 5.58106146679532777 - z + (z - 0.5) * log(z + 6.5)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#test_set = 1, 500, 120, 1800
#test_set = 10, 10, 20, 100
#test_set = 1, 5, 8, 15
#test_set = 2, 9, 10, 19
#test_set = 2, 12, 52, 962
#test_set = 79, 268, 195, 195+3723

#f = FisherExactTest()
#f.print_report(*test_set)
