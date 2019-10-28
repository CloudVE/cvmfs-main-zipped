/* dist_mat -- create a distance matrix in PHYLIP format for pairs of
*  specified individuals, including by default the reference sequence
*
*  argv[1] -- a Galaxy SNP table
*  argv[2] -- min coverage
*  argv[3] -- min quality
*  argv[4] -- name of reference species (or "none")
*  argv[5] -- 0 = distance from coverage; 1 = distance from genotype
*  argv[6] -- name of file for the numbers of informative SNPs
*  argv[7] -- name of file to write the Mega-format distance matrix
*  argv[k] for k > 7 have the form "13:fred", meaning that the 13th and 14th
*    columns (base 0) give the allele counts for the individual or group named
*    "fred".

What it does on Galaxy
This tool uses the selected SNP table to determine a "genetic distance" between each pair of selected individuals; the table of pairwise distances can be used by the Neighbor-Joining methods to construct a tree that depicts how the individuals are related. For a given pair of individuals, we find all SNP positions where both individuals have at least a minimum number of sequence "reads"; the individuals' distance at that SNP is defined as the absolute value of difference in the frequency of the first allele (equivalently: the second allele). For instance, if the first individuals has 5 reads of each allele and the second individual has respectivley 3 and 6 reads, then the frequencies are 1/2 and 1/3, giving a distance 1/6 at that SNP (provided that the minimum read total is at most 9). The output includes a report of the numbers of SNPs passing that thresold for each pair of individuals.

*/

#include "lib.h"

// bounds line length for a line of the Galaxy table

#define MOST 5000
#define MIN_SNPS 3

struct argument {
	int column;
	char *name;
} A[MOST];
int nA;	// number of individuals or groups + 1 (for the reference species)

#define MOST_INDIVIDUALS 100
#define SIZ 1+MOST_INDIVIDUALS // includes the reference

double tot_diff[SIZ][SIZ];
int ndiff[SIZ][SIZ], X[MOST];

int main(int argc, char **argv) {
	FILE *fp, *gp, *mega;
	char *p, *z = "\t\n", buf[MOST], name[100], B[100], C[100], D[100],
	  *nucs = "ACGT";
	int i, j, m, n, min_coverage, too_few, ref_allele = -1, has_ref,
	  min_quality, genotype;
	double fi, fj, dist;

	if (argc < 8)
		fatal("args: Galaxy-table min-cov min-qual min-snp ref-name genotype dist-out mega-out 13:fred 16:mary ...");
	min_coverage = atoi(argv[2]);
	min_quality = atoi(argv[3]);
	if (min_coverage <= 0 && min_quality <= 0)
		fatal("coverage and/or quality of SNPs should be constrained");

	if (same_string(argv[4], "none"))
		has_ref = 0;
	else {
		has_ref = 1;
		A[0].name = copy_string(argv[4]);
	}
	genotype = atoi(argv[5]);
	gp = ckopen(argv[6], "w");
	mega = ckopen(argv[7], "w");
	fprintf(mega, "#mega\n!Title: Galaxy;\n");
	  
	for (nA = has_ref, i = 8; i < argc; ++i, ++nA) {
		if (nA >= SIZ)
			fatal("Too many individuals");
		if (sscanf(argv[i], "%d:%s", &(A[nA].column), name) != 2)
			fatalf("bad arg: %s", argv[i]);
		A[nA].name = copy_string(name);
	}
	fprintf(mega,
	  "!Format DataType=Distance DataFormat=LowerLeft NTaxa=%d;\n\n",
	  nA);
	for (i = 0; i < nA; ++i)
		fprintf(mega, "[%d] #%s\n", i+1, A[i].name);
	fprintf(mega, "\n\n\n[");
	for (i = 1; i <= nA; ++i)
		fprintf(mega, "%4d", i);
	fprintf(mega, " ]\n");
	fp = ckopen(argv[1], "r");
	while (fgets(buf, MOST, fp)) {
		if (buf[0] == '#')
			continue;
		if (has_ref) {
			// get the reference allele
			if (sscanf(buf, "%*s %*s %s %s %*s %*s %*s %s", B, C, D)
			    != 3)
				fatalf("3 fields: %s", buf);
			if (strchr(nucs, B[0]) == NULL ||
			    strchr(nucs, C[0]) == NULL)
				fatalf("not nucs : %s %s", B, C);
			if (D[0] == B[0])
				ref_allele = 1;
			else if (D[0] == C[0])
				ref_allele = 2;
			else if (strchr(nucs, D[0]) != NULL)
				ref_allele = 3;
			else {
				if (D[0] != '-' && D[0] != 'N')
					fatalf("what is this: %s", D);
				ref_allele = -1;
			}
		}
			
		// X[i] = atoi(i-th word base-1)
		for (i = 1, p = strtok(buf, z); p != NULL;
		     ++i, p = strtok(NULL, z))
			X[i] = atoi(p);
		for (i = has_ref; i < nA; ++i) {
			m = A[i].column;
			if (X[m] + X[m+1] < min_coverage ||
			    X[m+3] < min_quality)
				continue;

			// frequency of the second allele
			if (genotype) {
				if (X[m+2] == -1)
					continue;	// no genotype
				fi = (double)X[m+2];
			} else
				fi = (double)X[m+1] / (double)(X[m]+X[m+1]);
			if (has_ref && ref_allele > 0) {
				ndiff[0][i]++;
				// reference allele might be different from both
				if (ref_allele == 1)
					tot_diff[0][i] += fi;
				else if (ref_allele == 2)
					tot_diff[0][i] += (1.0 - fi);
				else
					tot_diff[0][i] += 1.0;
			}
			for (j = i+1; j < nA; ++j) {
				n = A[j].column;
				if (X[n] + X[n+1] < min_coverage ||
				   X[n+3] < min_quality)
					continue;
				if (genotype && X[n+2] == -1)
					continue;
				ndiff[i][j]++;
				if (genotype)
					fj = (double)X[n+2];
				else
					fj = (double)X[n+1] /
					     (double)(X[n] + X[n+1]);
				fj -= fi;
				// add abs. value of difference in frequencies
				tot_diff[i][j] += (fj >= 0.0 ? fj : -fj);
			}

		}
	}
	for (i = too_few = 0; i < nA; ++i)
		for (j = i+1; j < nA; ++j)
			if (ndiff[i][j] < MIN_SNPS) {
				too_few = 1;
				fprintf(stderr,
				  "%s and %s have only %d informative SNPs\n",
				  A[i].name, A[j].name, ndiff[i][j]);
			}
	if (too_few)
		fatal("remove individuals or relax constraints");
		
	// print distances
	printf("%d\n", nA);
	for (i = 0; i < nA; ++i) {
		printf("%9s", A[i].name);
		fprintf(mega, "[%d] ", i+1);  
		for (j = 0; j < i; ++j) {
			dist = tot_diff[j][i]/(double)ndiff[j][i];
			printf(" %6.4f", dist);
			fprintf(mega, " %6.4f", dist);
		}
		fprintf(mega, "  \n");
		printf(" 0.0000");
		for (j = i+1; j < nA; ++j)
			printf(" %6.4f",
			  tot_diff[i][j]/(double)ndiff[i][j]);
		putchar('\n');
	}
	fprintf(mega, "\n\n\n\n\n");
	fclose(mega);

	// print numbers of SNPs
	for (i = 0; i < nA; ++i) {
		fprintf(gp, "%9s", A[i].name);
		for (j = 0; j < i; ++j)
			fprintf(gp, " %8d", ndiff[j][i]);
		fprintf(gp, "        0");
		for (j = i+1; j < nA; ++j)
			fprintf(gp," %8d", ndiff[i][j]);
		putc('\n', gp);
	}

	return 0;
}
