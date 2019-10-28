/* pop -- add four columns (allele counts, genotype, maximum quality) for a
*  specified population to a Galaxy SNP table, or enforce bounds
*
*  argv[1] = file containing a Galaxy table
*  argv[2] = lower bound on total coverage (-1 = no lower bound)
*  argv[3] = upper bound on total coverae (-1 if no bound)
*  argv[4] = lower bound on individual coverage (-1 = no bound)
*  argv[5] = lower bound on individual quality value (-1 = no bound)
*  argv[6] ... are the starting columns (base-1) for the chosen individuals

What it does on Galaxy
The user specifies that some of the individuals in the selected SNP table are form a "population" that has been previously defined using the Galaxy tool to select individuals from a SNP table. One option is for the program to append four columns to the table, giving the total counts for the two alleles, the "genotype" for the population and the maximum quality value, taken over all indivuals in the population. If all defined genotypes in the population are 2 (agree with the reference), the population's genotype is 2; similarly for 0; otherwise the genoype is 1 (unless all individuals have undefined genotype, in which case it is -1.  The other option is to remove rows from the table for which the total coverage for the population is either too low or too high, and/or if the individual coverage or quality value is too low.
*/

#include "lib.h"

// most characters allowed in a row of the table
#define MOST 50000

// column for the relevant individuals/groups
int col[MOST];
int nI;

int main(int argc, char **argv) {
	FILE *fp;
	char *p, *z = "\t\n", buf[MOST], trash[MOST];
	int X[MOST], m, i, A, B, G, Q, lo, hi, indiv, qual, g, q;

	if (argc < 3)
		fatalf("args: SNP-table low high col1 col2 ...");

	lo = atoi(argv[2]);
	hi = atoi(argv[3]);
	indiv = atoi(argv[4]);
	qual = atoi(argv[5]);
	for (i = 6, nI = 0; i < argc; ++i, ++nI)
		col[nI] = atoi(argv[i]);

	fp = ckopen(argv[1], "r");
	while (fgets(buf, MOST, fp)) {
		if (buf[0] == '#')
			continue;
		strcpy(trash, buf);
		// set X[i] = atoi(i-th word of s), i is base 0
		for (i = 1, p = strtok(trash, z); p != NULL;
		  ++i, p = strtok(NULL, z))
			X[i] = atoi(p);
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
