/* aggregate -- add four columns (allele counts, genotype, maximum quality) for
*  a specified population to a Galaxy SNP table
*
*  argv[1] = file containing a Galaxy table
*  argv[2] ... are the starting columns (base-1) for the chosen individuals

What it does on Galaxy
The user specifies that some of the individuals in a gd_snp dataset form a "population", by supplying a list that has been previously created using the Specify Individuals tool. The program appends a new "entity" (set of four columns) to the gd_snp table, analogous to the columns for an individual but containing summary data for the population as a group. These four columns give the total counts for the two alleles, the "genotype" for the population, and the maximum quality value, taken over all individuals in the population. If all defined genotypes in the population are 2 (agree with the reference), then the population's genotype is 2, and similarly for 0; otherwise the genotype is 1 (unless all individuals have undefined genotype, in which case it is -1).
*/

#include "lib.h"

// most characters allowed in a row of the table
#define MOST 5000

// column for the relevant individuals/groups
int col[MOST];
int nI;

int main(int argc, char **argv) {
	FILE *fp;
	char *p, *z = "\t\n", buf[MOST], trash[MOST];
	int X[MOST], m, i, A, B, G, Q, g;

	if (argc < 3)
		fatalf("args: SNP-table individual1 ...");

	for (i = 2, nI = 0; i < argc; ++i, ++nI)
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
			A += X[m];
			B += X[m+1];
			g = X[m+2];
			if (g != -1) {
				if (G == -1)	// first time
					G = g;
				else if (G != g)
					G = 1;
			}
			Q = MAX(Q, X[m+3]);
		}
		if (i < nI)	// check bounds on the population's individuals
			continue;
		// add columns
		if ((p = strchr(buf, '\n')) != NULL)
			*p = '\0';
		printf("%s\t%d\t%d\t%d\t%d\n", buf, A, B, G, Q);
	}

	return 0;
}
