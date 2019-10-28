/* filter_snps -- enforce constraints on a file of SNPs.
*
*  argv[1] = file containing a Galaxy table
*  argv[2] = 1 for a gd_snp file, 0 for gd_genotype
*  argv[3] = lower bound on total coverage (< 0 means interpret as percentage)
*  argv[4] = upper bound on total coveraae (< 0 means interpret as percentage;
*	     0 means no bound)
*  argv[5] = lower bound on individual coverage 
*  argv[6] = lower bound on individual quality value
*  argv[7] = lower bound on the number of defined genotypes
*  argv[8] = lower bound on the spacing between SNPs
*  argv[9] = reference-chromosome column (base-1); ref position in next column
*	If argv[8] == 0, then argv[7] must be 0.
*  argv[10] ... are the starting columns (base-1) for the chosen individuals.
*	If argc == 10, then only the lower bound on spacing between SNPs is
*	enforced.

What it does on Galaxy
For a gd_snp dataset, the user specifies that some of the individuals form a "population", by supplying a list that has been previously created using the Specify Individuals tool. SNPs are then discarded if their total coverage for the population is too low or too high, or if their coverage or quality score for any individual in the population is too low.

The upper and lower bounds on total population coverage can be specified either as read counts or as percentiles (e.g. "5%", with no decimal places). For percentile bounds the SNPs are ranked by read count, so for example, a lower bound of "10%" means that the least-covered 10% of the SNPs will be discarded, while an upper bound of, say, "80%" will discard all SNPs above the 80% mark, i.e. the top 20%. The threshold for the lower bound on individual coverage can only be specified as a plain read count.

For either a gd_snp or gd_genotype dataset, the user can specify a minimum number of defined genotypes (i.e., not -1) and/or a minimum spacing relative to the reference sequence. An error is reported if the user requests a minimum spacing but no reference sequence is available.

*/

#include "lib.h"

// most characters allowed in a row of the table
#define MOST 50000

// larger than any possible coverage
#define INFINITE_COVERAGE 100000000

char buf[MOST], chr[100], old_chr[100];

// column for the relevant individuals/groups
int col[MOST], *T;
int nI, lo, hi, min_space, min_geno, chr_col, pos, old_pos, gd_snp, X[MOST];

void get_X() {
	char *p, *z = "\t\n", trash[MOST];
	int i;

	strcpy(trash, buf);
	// set X[i] = atoi(i-th word of s), i is base 0
	for (i = 1, p = strtok(trash, z); p != NULL;
	  ++i, p = strtok(NULL, z)) {
		if (chr_col && i == chr_col)
			strcpy(chr, p);
		X[i] = atoi(p);
		if (chr_col && i == chr_col+1)
			pos = X[i];
	}
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

void OK() {
	if (chr_col == 0 || !same_string(chr, old_chr) ||
	    pos >= old_pos + min_space) {
		printf("%s", buf);
		if (chr_col) {
			strcpy(old_chr, chr);
			old_pos = pos;
		}
	}
}
int main(int argc, char **argv) {
	FILE *fp;
	int m, i, cov, tot_cov, indiv, qual, ngeno;

	if (argc < 10)
		fatalf("args: SNP-table gd_snp low-tot high-tot low-cov low-qual low-genotype low-space chr-col col1 col2 ...");

	for (i = 10, nI = 0; i < argc; ++i, ++nI)
		col[nI] = atoi(argv[i]);
	gd_snp = atoi(argv[2]);
	lo = atoi(argv[3]);
	hi = atoi(argv[4]);
	if (hi == 0)
		hi = INFINITE_COVERAGE;
	if (lo < 0 || hi < 0)
		find_lo(argv[1]);
	indiv = atoi(argv[5]);
	qual = atoi(argv[6]);
	min_geno = atoi(argv[7]);
	min_space = atoi(argv[8]);
	chr_col = atoi(argv[9]);

	// reality checks
	if (!gd_snp &&
	   (lo != 0 || hi != INFINITE_COVERAGE || indiv != 0 || qual != 0))
		fatal("cannot bound coverage or quality in gd_genotype file");
	if (chr_col == 0 && min_space != 0)
		fatalf("Specification of minimum spacing requires a reference sequence");
	if (indiv < 0 || qual < 0)
		fatalf("percentiles not implemented for individuals");
	
	// scan the SNPs
	fp = ckopen(argv[1], "r");
	while (fgets(buf, MOST, fp)) {
		if (buf[0] == '#')
			continue;
		get_X();
		for (i = tot_cov = ngeno = 0; i < nI; ++i) {
			m = col[i];
			if (gd_snp) {
				cov = (X[m]+X[m+1]);
				if (cov < indiv || X[m+3] < qual)
					break;
				tot_cov += cov;
			}
			if (X[m+2] != -1) 
				++ngeno;
		}
		if (i < nI)	// check bounds on the population's individuals
			continue;
		if (gd_snp && (tot_cov < lo || tot_cov > hi))
			continue;
		if (ngeno >= min_geno)
			// OK, except possibly for lower bound on spacing
			OK();
	}

	return 0;
}
