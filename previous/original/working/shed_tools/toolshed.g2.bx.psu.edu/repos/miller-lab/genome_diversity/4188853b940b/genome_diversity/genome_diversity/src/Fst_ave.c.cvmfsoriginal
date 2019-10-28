/* Fst_ave -- determine four FST values between two specified populations,
*  and optionally between several pairs of random populations
*
*    argv{1] = a Galaxy SNP table. For each of several individuals, the table
*              has four columns (#A, #B, genotype, quality).
*    argv[2] = 1 if FST is estimated from SAMtools genotypes; 0 means use
*	        read-coverage data.
*    argv[3] = lower bound, for individual quality value if argv[2] = 1,
*	       or for total number of reads per population if argv[2] = 0.
*	       SNPs not satisfying these lower bounds are ignored.
*    argv[4] = 1 to discard SNPs that appear fixed in the two populations
*    argv[5] = k says report the maximum and average FST over k randomly
*              chosen splits into two populations of two original sizes
*    argv[6], argv[7], ...,  have the form "13:1", "13:2" or "13:0", meaning
*             that the 13th and 14th columns (base 1) give the allele counts
*             (and column 15 gives the genotype) for an individual that is in
*	       population 1, in population 2, or in neither population.

What it does on Galaxy

The user specifies a SNP table and two "populations" of individuals, both previously defined using the Galaxy tool to specify individuals from a SNP table. No individual can be in both populations. Other choices are as follows.

Data source. The allele frequencies of a SNP in the two populations can be estimated either by the total number of reads of each allele, or by adding the frequencies inferred from genotypes of individuals in the populations.

After specifying the data source, the user sets lower bounds on amount of data required at a SNP. For estimating the FST using read counts, the bound is the minimum count of reads of the two alleles in a population. For estimations based on genotype, the bound is the minimum reported genotype quality per individual. SMPs not meeting these lower bounds are ignored.

The user specifies whether SNPs where both populations appear to be fixed for the same allele should be retained or discarded.

Finally, the user decides whether to use randomizations. If so, then the user specifies how many randomly generated population pairs (retaining the numbers of individuals of the originals) to generate, as well as the "population" of additional individuals (not in the first two populations) that can be used in the randomization process.

The program prints the following measures of FST for the two populations.
1. The formulation by Sewall Wright (average over FSTs for all SNPs).
2. The Weir-Cockerham estimator (average over FSTs for all SNPs).
3. The Reich-Patterson estimator (average over FSTs for all SNPs).
4. The population-based Reich-Patterson estimator.

If randomizations were requested, it prints a summary for each of the four definitions of FST that includes the maximum and average value, and the highest-scoring population pair (if any scored higher than the two user-specified populations).

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

// maximum length of a line from the table
#define MOST 50000

// information about the specified individuals
// x is an array of nI values 0, 1, or 2;
// shuffling x creates random "populations"
int col[MOST], x[MOST];
int nI, lower_bound, discard, genotypes, nsnp, nfail;
double F_wright, F_weir, F_reich, N_reich, D_reich;

// each SNP has an array of counts
struct count {
	int A, B;
};

// linked list summarizes the Galaxy table
struct snp {
	struct count *c;
	struct snp *next;
} *start, *last;

/* For each of wright, weir and reich, we observe allele counts A1 and A2
*  for one allele in the two populations, and B1 and B2 for the other allele.
*/

// given the two populations specified by x[], compute four corresponding FSTs
void pop_Fst() {
	double N, D;
	struct snp *s;
	int i, A1, B1, A2, B2, too_few;


	// scan the SNPs
	F_wright = F_weir = F_reich = N_reich = D_reich = 0.0;
	nsnp = nfail = 0;
	for (s = start; s != NULL; s = s->next) {
		// get counts for the two populations at this SNP
		for (A1 = B1 = A2 = B2 = i = 0; i < nI; ++i) {
			if (s->c[i].A < 0) // no genotypes
				continue;
			if (x[i] == 1) {
				A1 += s->c[i].A;
				B1 += s->c[i].B;
			} else if (x[i] == 2) {
				A2 += s->c[i].A;
				B2 += s->c[i].B;
			}
		}
		if (discard && ((A1 == 0 && A2 == 0) || (B1 == 0 && B2 == 0)))
			continue;	// fixed in these two populations
		too_few = (genotypes ? 1 : lower_bound);
		if (A1+B1 >= too_few && A2+B2 >= too_few) {
			++nsnp;
			wright(A1, A2, B1, B2, &N, &D);
			if (D != 0.0)
				F_wright += N/D;
			else
				++nfail;
			weir(A1, A2, B1, B2, &N, &D);
			if (D != 0.0)
				F_weir += N/D;
			else
				++nfail;
			reich(A1, A2, B1, B2, &N, &D);
			N_reich += N;
			D_reich += D;
			if (D != 0.0)
				F_reich += N/D;
			else
				++nfail;
		}
	}
	F_wright /= nsnp;
	F_weir /= nsnp;
	N_reich /= nsnp;
	D_reich /= nsnp;
	F_reich /= nsnp;
}

/* shuffle the values x[0], x[1], ... , x[nI-1];
*  Uses Algorithm P in page 125 of "The Art of Computer Programming (Vol II)
*  Seminumerical Programming", by Donald Knuth, Addison-Wesley, 1971.
*/
void shuffle() {
	int i, j, temp;

	for (i = nI - 1; i > 0; --i) {
		// swap what's in location i with location j, where 0 <= j <= i
		j = random() % (i+1);
		temp = x[i];
		x[i] = x[j];
		x[j] = temp;
	} 
}

int main(int argc, char **argv) {
	FILE *fp;
	char *p, *z = "\t\n", buf[MOST];
	int X[MOST], nshuff, n, i, j, k, saw[3], larger0, larger1, larger2,
	  larger3, best_x0[MOST], best_x1[MOST], best_x2[MOST], best_x3[MOST];
	struct snp *new;
	double F, F0, F1, F2, F3, tot_F0, tot_F1, tot_F2, tot_F3,
	  largest_F0, largest_F1, largest_F2, largest_F3;

	if (argc < 7)
		fatal("args: table data-source lower_bound discard? #shuffles n:1 m:2 ...");

	// handle command-line arguments
	genotypes = atoi(argv[2]);
	lower_bound = atoi(argv[3]);
	if (!genotypes && lower_bound <= 0)
		fatal("minimum coverage should exceed 0");
	discard = atoi(argv[4]);
	nshuff = atoi(argv[5]);
	saw[0] = saw[1] = saw[2] = 0;
	// populations 1 and 2 must be disjoint
	for (i = 6; i < argc; ++i) {
		if (sscanf(argv[i], "%d:%d", &j, &k) != 2)
			fatalf("not like 13:2 : %s", argv[i]);
		if (k < 0 || k > 2)
			fatalf("not population 0, 1 or 2: %s", argv[i]);
		saw[k] = 1;
		// seen this individual (i.e., column) before??
		for (n = 0; n < nI && col[n] != j; ++n)
			;
		if (n < nI) { // OK if one of the populations is 0
			if (k > 0) {
				if (x[n] > 0 && x[n] != k)
				  fatalf("column %d is in both populations", j);
				x[n] = k;
			}
		} else {
			col[nI] = j;
			x[nI] = k;
			++nI;
		}
	}
	if (saw[1] == 0)
		fatal("population 1 is empty");
	if (saw[2] == 0)
		fatal("population 2 is empty");

	// read the table of SNPs and store the essential allele counts
	fp = ckopen(argv[1], "r");
	while (fgets(buf, MOST, fp)) {
		if (buf[0] == '#')
			continue;
		new = ckalloc(sizeof(*new));
		new->next = NULL;
		new->c = ckalloc(nI*sizeof(struct count));
		// set X[i] = atoi(i-th word of buf), i is base 1
		for (i = 1, p = strtok(buf, z); p != NULL;
		  ++i, p = strtok(NULL, z))
			X[i] = atoi(p);
		for (i = 0; i < nI; ++i) {
			n = col[i];
			if (genotypes) {
				k = X[n+2];
				if (k == -1)
					new->c[i].A = new->c[i].B = -1;
				else {
					new->c[i].A = k;
					new->c[i].B = 2 - k;
				}
			} else {
				new->c[i].A = X[n];
				new->c[i].B = X[n+1];
			}
		}
		if (start == NULL)
			start = new;
		else
			last->next = new;
		last = new;
	}
	fclose(fp);

	pop_Fst();
	printf("Using %d SNPs, we compute:\n", nsnp);
	printf("Average Reich-Patterson FST is %5.5f.\n", F2 = F_reich);
	printf("The population-based Reich-Patterson Fst is %5.5f.\n",
	  F3 = N_reich/D_reich);
	printf("Average Weir-Cockerham FST is %5.5f.\n", F1 = F_weir);
	printf("Average Wright FST is %5.5f.\n", F0 = F_wright);
	if (nfail > 0)
		printf("WARNING: %d of %d FSTs could not be computed\n",
		  nfail, 3*nsnp);
	if (nshuff == 0)
		return 0;

	// do the following only if randomization is requested
	for (j = 0; j < nI; ++j)
		best_x0[j] = best_x1[j] = best_x2[j] = best_x3[j] = x[j];
	tot_F0 = tot_F1 = tot_F2 = tot_F3 =
	  largest_F0 = largest_F1 = largest_F2 = largest_F3 = 0.0;
	larger0 = larger1 = larger2 = larger3 = 0;
	for (i = 0; i < nshuff; ++i) {
		shuffle();
		pop_Fst();

		// Wright
		if ((F = F_wright) > F0)
			++larger0;
		if (F > largest_F0) {
			largest_F0 = F;
			for (j = 0; j < nI; ++j)
				best_x0[j] = x[j];
		}
		tot_F0 += F;
/*
		if (all)	// make this optional?
			printf("%d: %f\n", i+1, F);
*/

		// Weir
		if ((F = F_weir) > F1)
			++larger1;
		if (F > largest_F1) {
			largest_F1 = F;
			for (j = 0; j < nI; ++j)
				best_x1[j] = x[j];
		}
		tot_F1 += F;

		// Riech average
		if ((F = F_reich) > F2)
			++larger2;
		if (F > largest_F2) {
			largest_F2 = F;
			for (j = 0; j < nI; ++j)
				best_x2[j] = x[j];
		}
		tot_F2 += F;

		// Reich population
		if ((F = (N_reich/D_reich)) > F3)
			++larger3;
		if (F > largest_F3) {
			largest_F3 = F;
			for (j = 0; j < nI; ++j)
				best_x3[j] = x[j];
		}
		tot_F3 += F;
	}
	printf("\nOf %d random groupings:\n", nshuff);
	printf("%d had a larger average Wright FST (max %5.5f, mean %5.5f)\n",
	  larger0, largest_F0, tot_F0/nshuff);
	if (largest_F0 > F0) {
		printf("first columns for the best two populations:\n");
		for (i = 0; i < nI; ++i)
			if (best_x0[i] == 1)
				printf("%d ", col[i]);
		printf("and\n");
		for (i = 0; i < nI; ++i)
			if (best_x0[i] == 2)
				printf("%d ", col[i]);
		putchar('\n');
		putchar('\n');
	}
	printf("%d had a larger average Weir-Cockerham FST (max %5.5f, mean %5.5f)\n",
	  larger1, largest_F1, tot_F1/nshuff);
	if (largest_F1 > F1) {
		printf("first columns for the best two populations:\n");
		for (i = 0; i < nI; ++i)
			if (best_x1[i] == 1)
				printf("%d ", col[i]);
		printf("and\n");
		for (i = 0; i < nI; ++i)
			if (best_x1[i] == 2)
				printf("%d ", col[i]);
		putchar('\n');
		putchar('\n');
	}
	printf("%d had a larger average Reich-Patterson FST (max %5.5f, mean %5.5f)\n",
	  larger2, largest_F2, tot_F2/nshuff);
	if (largest_F2 > F2) {
		printf("first columns for the best two populations:\n");
		for (i = 0; i < nI; ++i)
			if (best_x2[i] == 1)
				printf("%d ", col[i]);
		printf("and\n");
		for (i = 0; i < nI; ++i)
			if (best_x2[i] == 2)
				printf("%d ", col[i]);
		putchar('\n');
		putchar('\n');
	}
	printf("%d had a larger Reich-Patterson population FST (max %5.5f, mean %5.5f)\n",
	  larger3, largest_F3, tot_F3/nshuff);
	if (largest_F3 > F3) {
		printf("first columns for the best two populations:\n");
		for (i = 0; i < nI; ++i)
			if (best_x3[i] == 1)
				printf("%d ", col[i]);
		printf("and\n");
		for (i = 0; i < nI; ++i)
			if (best_x3[i] == 2)
				printf("%d ", col[i]);
		putchar('\n');
		putchar('\n');
	}

	return 0;
}
