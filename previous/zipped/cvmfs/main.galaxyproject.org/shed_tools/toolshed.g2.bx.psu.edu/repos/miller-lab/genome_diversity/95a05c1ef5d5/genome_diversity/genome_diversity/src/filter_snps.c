/* pop -- add four columns (allele counts, genotype, maximum quality) for a
*  specified population to a Galaxy SNP table, or enforce bounds
*
*  argv[1] = file containing a Galaxy table
*  argv[2] = lower bound on total coverage (> 0 means interpret as percentage)
*  argv[3] = upper bound on total coverae (> 0 means interpret as percentage)
*  argv[4] = lower bound on individual coverage 
*  argv[5] = lower bound on individual quality value
*  argv[6] ... are the starting columns (base-1) for the chosen individuals

What it does on Galaxy
The user specifies that some of the individuals in a gd_snp dataset form a "population", by supplying a list that has been previously created using the Specify Individuals tool. SNPs are then discarded if their total coverage for the population is too low or too high, or if their coverage or quality score for any individual in the population is too low.

*/

#include "lib.h"

// most characters allowed in a row of the table
#define MOST 5000

char buf[MOST];

// column for the relevant individuals/groups
int col[MOST], *T;
int nI, lo, hi, X[MOST];

void get_X() {
	char *p, *z = "\t\n", trash[MOST];
	int i;

	strcpy(trash, buf);
	// set X[i] = atoi(i-th word of s), i is base 0
	for (i = 1, p = strtok(trash, z); p != NULL;
	  ++i, p = strtok(NULL, z))
		X[i] = atoi(p);
}

int compar(int *a, int *b) {
	if (*a < *b)
		return -1;
	if (*a > *b)
		return 1;
	return 0;
}

void find_lo(char *filename) {
	FILE *fp = ckopen(filename, "r");
	int n, m, i, k;

	for (n = 0; fgets(buf, MOST, fp); ++n)
		;
	T = ckalloc(n*sizeof(int));
	rewind(fp);
	for (k = 0; fgets(buf, MOST, fp); ++k) {
		get_X();
		for (i = T[k] = 0; i < nI; ++i) {
			m = col[i];
			T[k] += (X[m]+X[m+1]);
		}
	}
	qsort((void *)T, (size_t)n, sizeof(int), (const void *)compar);
	if (lo < 0) {
		lo = -lo;
		if (lo > 100)
			fatal("cannot have low > 100%");
		lo = T[(n*lo)/100];
	}
	if (hi < 0) {
		hi = -hi;
		if (hi > 100)
			fatal("cannot have high > 100%");
		hi = T[(n*hi)/100];
	}
	free(T);
	fclose(fp);
}

int main(int argc, char **argv) {
	FILE *fp;
	char *p;
	int m, i, A, B, G, Q, indiv, qual, g, q;

	if (argc < 3)
		fatalf("args: SNP-table low high col1 col2 ...");

	for (i = 6, nI = 0; i < argc; ++i, ++nI)
		col[nI] = atoi(argv[i]);
	lo = atoi(argv[2]);
	hi = atoi(argv[3]);
	if (hi == 0)
		hi = 100000000;
	if (lo < 0 || hi < 0)
		find_lo(argv[1]);
	indiv = atoi(argv[4]);
	qual = atoi(argv[5]);
	if (indiv < 0 || qual < 0)
		fatalf("percentiles not implemented for individuals");

	fp = ckopen(argv[1], "r");
	while (fgets(buf, MOST, fp)) {
		if (buf[0] == '#')
			continue;
		get_X();
		for (i = A = B = Q = 0, G = -1; i < nI; ++i) {
			m = col[i];
			if (X[m]+X[m+1] < indiv || (q = X[m+3]) < qual)
				break;
			A += X[m];
			B += X[m+1];
			g = X[m+2];
			if (g != -1) {
				if (G == -1)	// first time
					G = g;
				else if (G != g)
					G = 1;
			}
			Q = MAX(Q, q);
		}
		if (i < nI)	// check bounds on the population's individuals
			continue;
		if (lo == -1 && hi == -1 && indiv == -1 && qual == -1) {
			// add columns
			if ((p = strchr(buf, '\n')) != NULL)
				*p = '\0';
			printf("%s\t%d\t%d\t%d\t%d\n", buf, A, B, G, Q);
		} else if (A+B >= lo && (hi == -1 || A+B <= hi))
			// coverage meets the population-level restrictions
			printf("%s", buf);
	}

	return 0;
}
