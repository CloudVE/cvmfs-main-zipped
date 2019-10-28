/* admix_prep -- prepare the ".ped" and ".map" files (PLINK format) for input to
*  the "admixture" program.
*
*  argv[1] -- a Galaxy SNP table
*  argv[2] -- required number of reads for each individual to use a SNP
*  argv[3] -- required genotype quality for each individual to use a SNP
*  argv[4] -- minimum spacing between SNPs on the same scaffold
*  argv[k] for k > 4 have the form "13:fred", meaning that the 13th and 14th
*    columns (base 0) give the allele counts for the individual or group named
*    "fred".

What it does on Galaxy
The tool converts a SNP table into two tables, called "admix.map" and "admix.ped", needed for estimating the population structure. The user can read or download those files, or simply pass this tool's output on to other programs. The user imposes conditions on which SNPs to consider, such as the minimum coverage and/or quality value for every individual, or the distance to the closest SNP in the same contig (as named in the first column of the SNP table). A useful piece of information produced by the tool is the number of SNPs meeting those conditions, which can be found by clicking on the "eye" after the program runs.

*/

#include "lib.h"

// bounds line length for a line of the Galaxy table
#define MOST 50000
struct individual {
	int column;
	char *name;
} I[MOST/8]; // each individual has 4 columns and 4 tab characters
int nI;	// number of individuals
int X[MOST];	// integer values in a row of the SNP table

// bounds the number of SNPs that can be kept
#define MAX_KEEP 10000000
char *S[MAX_KEEP];	// S[i] is a row of 2*nI alleles
int nK;

int main(int argc, char **argv) {
	FILE *fp, *ped, *map;
	char *p, *z = " \t\n", buf[MOST], trash[MOST], name[100], *s,
	  scaf[100], prev_scaf[100];
	int i, j, m, min_coverage, min_quality, min_space, nsnp, genotype,
	   pos, prev_pos;

	if (argc < 5)
		fatal("args: Galaxy-table min-cov min-qual min-space 13:fred 16:mary ...");
	min_coverage = atoi(argv[2]);
	min_quality = atoi(argv[3]);
	min_space = atoi(argv[4]);

	for (i = 5; i < argc; ++i, ++nI) {
		if (nI >= MOST/8)
			fatal("Too many individuals");
		if (sscanf(argv[i], "%d:%s", &(I[nI].column), name) != 2)
			fatalf("bad arg: %s", argv[i]);
		I[nI].name = copy_string(name);
	}

	map = ckopen("admix.map", "w");

	fp = ckopen(argv[1], "r");
	prev_scaf[0] = '\0';
	prev_pos = 0;
	for (nsnp = 0; fgets(buf, MOST, fp); ) {
		if (buf[0] == '#')
			continue;
		++nsnp;
		if (sscanf(buf, "%s %d", scaf, &pos) != 2)
			fatalf("choke: %s", buf);
		if (same_string(scaf, prev_scaf)) {
			if (pos < prev_pos + min_space)
				continue;
		} else {
			strcpy(prev_scaf, scaf);
			prev_pos = -min_space;
		}

		// X[i] = atoi(i-th word base-1)
		strcpy(trash, buf);
		for (i = 1, p = strtok(trash, z); p != NULL;
		     ++i, p = strtok(NULL, z))
			X[i] = atoi(p);
		for (i = 0; i < nI; ++i) {
			m = I[i].column;
			if (X[m] + X[m+1] < min_coverage || X[m+3] < min_quality)
				break;
		}
		if (i < nI)
			continue;
		prev_pos = pos;
		
		if (nK >= MAX_KEEP)
			fatal("Too many SNPs");
		fprintf(map, "1 snp%d 0 %d\n", nsnp, nsnp+1);
		s = S[nK++] = ckalloc(2*nI*sizeof(char));
		for (i = j = 0; i < nI; ++i, j += 2) {
			genotype = X[I[i].column+2];
			if (genotype == 2)
				s[j] = s[j+1] = '1';
			else if (genotype == 0)
				s[j] = s[j+1] = '2';
			else if (genotype == 1) {
				s[j] = '1';
				s[j+1] = '2';
			} else	// undefined genotype
				s[j] = s[j+1] = '0';
		}
	}

	fclose(map);

	ped = ckopen("admix.ped", "w");
	for (i = 0; i < nI; ++i) {
		fprintf(ped, "%s 1 0 0 1 1", I[i].name);
		for (j = 0; j < nK; ++j)
			fprintf(ped, " %c %c", S[j][2*i], S[j][2*i+1]);
		putc('\n', ped);
	}

	printf("Using %d of %d SNPs\n", nK, nsnp);
	fclose(ped);

	return 0;
}
