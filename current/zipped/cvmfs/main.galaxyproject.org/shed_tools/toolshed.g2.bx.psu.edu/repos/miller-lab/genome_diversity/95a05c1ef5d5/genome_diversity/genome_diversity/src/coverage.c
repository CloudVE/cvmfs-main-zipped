/* coverage -- report distributions of SNP coverage or quality for individuals,
*  or coverage for populations
*
*    argv{1] -- a Galaxy SNP table. For each individuals, the table has four
*		columns (count of each allele, genotype, quality).
*    argv[2] -- 0 = sequence coverage, 1 = genotype quality
*    argv[3] -- file name for the text version of output (input for producing
*		the graphical summary goes to stdout)
*    argv[4], argv[5], ...,  have the form "13:fred",  meaning that the 13th
*		14th, and 16th columns (base 1) give the two allele counts
*		and the quality for "fred", where "fred" can be the name of
*		a population with several individuals (all named "fred")
What it does on Galaxy
The tool reports distributions of SNP reliability indicators for individuals or populations. The reliability can be measured by either the sequence coverage or the SAMtools quality value, though the notion of a population-level quality is not supported. Textual and graphical reports are generated, where the text output gives the cumulative distributions.
*/

#include "lib.h"

// maximum length of a line from the table
#define MOST 5000

// the largest coverage or quality value being considered
#define MAX_VAL 1000

FILE *gp;	// for text output

// a population is the set of all indivuals with the same name
// (perhaps just a single individual)
struct pop {
	int cov, n[MAX_VAL+1];
	long long sum, tot;
	char *name;
} P[MOST/4];
int nP;	// number of populations

// maps column to population
struct individual {
	int col, pop;
} I[MOST/4];
int nI;

/* Report the distribution for each individual. P[i].n[k] is the number of SNPs
*  of value (coverage or quality) k in population i, for k < MAX_VAL;
*  I[i].n[MAX_VAL] is the number of SNPs of value k >= MAX_VAL.
*  We print the percentages, p, of SNPs with value <= k, ending when all
*  populations have reached a p >= 98%.
*/
void print_cov() {
	int i, j, k, last_j;
	long long sum;

	// find where to stop printing
	for (last_j = i = 0; i < nP; ++i) {
		for (sum = j = 0; j <= MAX_VAL; ++j)
			sum += P[i].n[j];
		P[i].tot = sum;
		for (sum = j = 0; j <= MAX_VAL; ++j) {
			sum += P[i].n[j];
			if (sum >= 0.98*P[i].tot)
				break;
		}
		last_j = MAX(last_j, j);
	}


	++last_j;
	// print to stdout the output for graphing; not broken into short lines
	for (j = 0; j < last_j; ++j)
		printf("\t%3d", j);
	putchar('\n');
	for (i = 0; i < nP; ++i) {
		printf("%s", P[i].name);
		for (sum = j = 0; j < last_j; ++j) {
			sum += P[i].n[j];
			printf("\t%4.2f", 100.0*(float)sum/(float)P[i].tot);
		}
		putchar('\n');
	}

	// print a user-friendly version to the named file
	// <= 20 numbers per row
	for (j = 0; j < last_j; j += 20) {
		fprintf(gp, "\n          ");
		for (k = j; k < MIN(j+20, last_j); ++k)
			fprintf(gp, "%3d", k);
		for (i = 0; i < nP; ++i) {
			fprintf(gp, "\n%10s", P[i].name);
			for (k = j; k < MIN(j+20, last_j); ++k) {
				P[i].sum += P[i].n[k];
				fprintf(gp, "%3lld",
				  MIN(99, 100*P[i].sum/P[i].tot));
			}
		}
		fprintf(gp,"\n\n");
	}
}

int main(int argc, char **argv) {
	FILE *fp;
	char buf[MOST], *z = " \t\n", *p;
	int X[MOST], i, j, cov, m, quality, is_pop;

	if (argc < 5)
		fatal("args: SNP-file quality-value? out-name 13:fred ... ");
	quality = atoi(argv[2]);
	gp = ckopen(argv[3], "w");
	// record the individuals and populations
	for (nI = 0, i = 4; i < argc; ++i, ++nI) {
		if (nI >= MOST)
			fatal("Too many individuals");
		// allow spaces in names
		if ((p = strchr(argv[i], ':')) == NULL)
			fatalf("no colon: %s", argv[i]);
		I[nI].col = atoi(argv[i]);
		for (j = 0; j < nP && !same_string(p+1, P[j].name); ++j)
			;
		if (j == nP) { // new population
			is_pop = 1;
			P[nP++].name = copy_string(p+1);
		}
		I[nI].pop = j;
	}
	if (is_pop && quality)
		fatal("quality values for a population are not supported.");

	// Record the number of SNPs with coverage 0, 1, ..., MAX_VAL-1,
	// or >= MAX_VAL for each individual.
	fp = ckopen(argv[1], "r");
	while (fgets(buf, MOST, fp)) {
		if (buf[0] == '#')
			continue;
		// P[i].cov is the total coverage for all individuals in pop i
		for (i = 0; i < nP; ++i)
			P[i].cov = 0;
		// X[i] = atoi(i-th word base-1)
		for (i = 1, p = strtok(buf, z); p != NULL;
		     ++i, p = strtok(NULL, z))
			X[i] = atoi(p);
		for (i = 0; i < nI; ++i) {
			m = I[i].col;
			if (quality)
				cov = X[m+3];
			else
				cov = X[m] + X[m+1];
			P[I[i].pop].cov += cov;
		}
		for (i = 0; i < nP; ++i)
			P[i].n[MIN(P[i].cov, MAX_VAL)]++;
	}

	// Print the distributions.
	print_cov();

	return 0;
}
