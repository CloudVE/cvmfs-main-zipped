/* aggregate -- add four columns (allele counts, genotype, maximum quality) for
*  a specified population to a Galaxy SNP table
*
*  argv[1] = file containing a Galaxy table
*  argv[2] = 0 for a gd_genotype file, 1 for a gd_snp file
*  argv[3] ... are the starting columns (base-1) for the chosen individuals

What it does on Galaxy
The user specifies that some of the individuals in a gd_snp or gd_genotype dataset form a "population", by supplying a list that has been previously created using the Specify Individuals tool. The program appends a new "entity" (set of four columns for a gd_snp table, or one column for a gd_genotype table), analogous to the column(s) for an individual but containing summary data for the population as a group. For a gd_snp table, these four columns give the total counts for the two alleles, the "genotype" for the population, and the maximum quality value, taken over all individuals in the population. If all defined genotypes in the population are 2 (agree with the reference), then the population's genotype is 2, and similarly for 0; otherwise the genotype is 1 (unless all individuals have undefined genotype, in which case it is -1). For a gd_genotype file, only the aggregate genotype is appended.
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
	int X[MOST], m, i, A, B, G, Q, g, gd_snp;

	if (argc < 3)
		fatalf("args: SNP-table typedef individual1 ...");

	gd_snp = atoi(argv[2]);
	for (i = 3, nI = 0; i < argc; ++i, ++nI)
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
			if (gd_snp) {
				A += X[m];
				B += X[m+1];
				Q = MAX(Q, X[m+3]);
			}
			g = X[m+2];
			if (g != -1) {
				if (G == -1)	// first time
					G = g;
				else if (G != g)
					G = 1;
			}
		}
		if (i < nI)	// check bounds on the population's individuals
			continue;
		// add columns
		if ((p = strchr(buf, '\n')) != NULL)
			*p = '\0';
		if (gd_snp)
			printf("%s\t%d\t%d\t%d\t%d\n", buf, A, B, G, Q);
		else
			printf("%s\t%d\n", buf, G);
	}

	return 0;
}
