/* Fst_column -- add an Fst column to a Galaxy table
*
*    argv{1] = a Galaxy SNP table. For each of several individuals, the table
*              has four columns (#A, #B, genotype, quality).
*    argv[2] = 1 if Fst is estimated from SAMtools genotypes; 0 means use
*	        read-coverage data.
*    argv[3] = lower bound for total number of reads per population
*    argv[4] = lower bound for individual quality value
*    argv[5] = 1 to retain SNPs that fail to satisfy the lower bound and set
*	       Fst = -1; delete them if argv[4] = 0.
*    argv[6] = 1 to discard SNPs that appear fixed in the two populations
*    argv[7] = 0 for the original Wright form, 1 for Weir, 2 for Reich
*    argv[8], argv[9], ...,  have the form "13:1" or "13:2", meaning that
*             the 13th, 14th, and 15th columns (base 1) give the allele counts
*             and genotype for an individual that is in population 1 or
*	      population 2, respectively.

What It Does on Galaxy

The user specifies a SNP table and two "populations" of individuals, both previously defined using the Galaxy tool to specify individuals from a SNP table. No individual can be in both populations. Other choices are as follows.

Data source. The allele frequencies of a SNP in the two populations can be estimated either by the total number of reads of each allele, or by adding the frequencies inferred from genotypes of individuals in the populations.

After specifying the data source, the user sets lower bounds on amount of data required at a SNP. For estimating the Fst using read counts, the bound is the minimum count of reads of the two alleles in a population. For estimations based on genotype, the bound is the minimum reported genotype quality per individual.

The user specifies whether the SNPs that violate the lower bound should be ignored or the Fst set to -1.

The user specifies whether SNPs where both populations appear to be fixed for the same allele should be retained or discarded.

Finally, the user chooses which definition of Fst to use: Wright's original definition, the Weir-Cockerham unbiased estimator, or the Reich-Patterson estimator.

A column is appended to the SNP table giving the Fst for each retained SNP.

References:

Sewall Wright (1951) The genetical structure of populations. Ann Eugen 15:323-354.

B. S. Weir and C. Clark Cockerham (1984) Estimating F-statistics for the analysis of population structure. Evolution 38:1358-1370.

Weir, B.S. 1996. Population substructure. Genetic data analysis II, pp. 161-173. Sinauer Associates, Sunderland, MA.

David Reich, Kumarasamy Thangaraj, Nick Patterson, Alkes L. Price, and Lalji Singh (2009) Reconstructing Indian population history. Nature 461:489-494, especially Supplement 2.  

Their effectiveness for computing FSTs when there are many SNPs but few individuals is discussed in the followoing paper.

Eva-Maria Willing, Christine Dreyer, Cock van Oosterhout (2012) Estimates of genetic differentiation measured by FST do not necessarily require large sample sizes when using many SNP markers. PLoS One 7:e42649.

*/

#include "lib.h"
#include "Fst_lib.h"

// most characters allowed in a row of the table
#define MOST 50000

// column and population for the relevant individuals/groups
int col[MOST], pop[MOST];
int nI;

int main(int argc, char **argv) {
	FILE *fp;
	char *p, *z = "\t\n", buf[MOST], trash[MOST];
	int X[MOST], min_cov, min_qual, retain, discard, unbiased, genotypes,
	  n, i, g, A1, B1, A2, B2, saw[3], x1, y1, x2, y2;
	double F, N, D;

	if (argc < 7)
		fatal("args: table data-source lower-bound retain? discard? unbiased? n:1 m:2 ...");
	genotypes = atoi(argv[2]);
	min_cov = atoi(argv[3]);
	min_qual = atoi(argv[4]);
	retain = atoi(argv[5]);
	discard = atoi(argv[6]);
	unbiased = atoi(argv[7]);
	saw[1] = saw[2] = 0;
	for (i = 8; i < argc; ++i, ++nI) {
		if (sscanf(argv[i], "%d:%d", &(col[nI]), &(pop[nI])) != 2)
			fatalf("not like 13:2 : %s", argv[i]);
		if (pop[nI] < 1 || pop[nI] > 2)
			fatalf("not population 1 or 2: %s", argv[i]);
		saw[pop[nI]] = 1;
		// seen this individual before?
		for (n = 0; n < nI && col[n] != col[nI]; ++n)
			;
		if (n < nI)
			fatalf("individual at column %d is mentioned twice",
			  col[n]);
	}
	if (saw[1] == 0)
		fatal("population 1 is empty");
	if (saw[2] == 0)
		fatal("population 2 is empty");

	fp = ckopen(argv[1], "r");
	while (fgets(buf, MOST, fp)) {
		if (buf[0] == '#')
			continue;
		strcpy(trash, buf);
		// set X[i] = atoi(i-th word of s), i is base 0
		for (i = 1, p = strtok(trash, z); p != NULL;
		  ++i, p = strtok(NULL, z))
			X[i] = atoi(p);
		for (i = A1 = B1 = A2 = B2 = x1 = y1 = x2 = y2 = 0;
		     i < nI; ++i) {
			n = col[i];
			g = X[n+2];	// save genotype

			if (genotypes) {
				if (g == -1)
					continue;
			} else if (X[n+3] < min_qual)
				continue;
			if (pop[i] == 1) {
				// column n (base 1) corresponds to entry X[n]
				x1 += X[n];
				y1 += X[n+1];
				if (genotypes) {
					A1 += g;
					B1 += (2 - g);
				} else {
					A1 += X[n];
					B1 += X[n+1];
				}
			} else if (pop[i] == 2) {
				x2 += X[n];
				y2 += X[n+1];
				if (genotypes) {
					A2 += g;
					B2 += (2 - g);
				} else {
					A2 += X[n];
					B2 += X[n+1];
				}
			}
		}
		if (discard && ((A1 == 0 && A2 == 0) || (B1 == 0 && B2 == 0)))
			continue; // not variable in the two populations
		if (!genotypes && (x1+y1 < min_cov || x2+y2 < min_cov))
			F = -1.0;
		else {
			if (unbiased == 0)
				wright(A1, A2, B1, B2, &N, &D);
			else if (unbiased == 1)
				weir(A1, A2, B1, B2, &N, &D);
			else if (unbiased == 2)
				reich(A1, A2, B1, B2, &N, &D);
			else
				fatal("impossible value of 'unbiased'");
			if (D == 0.0)
				continue;	// ignore these SNPs
			else
				F = N/D;
		}
		if (F == -1.0 && !retain)
			continue;
		if ((p = strchr(buf, '\n')) != NULL)
			*p = '\0';
		printf("%s\t%5.4f\n", buf, F);
	}

	return 0;
}
