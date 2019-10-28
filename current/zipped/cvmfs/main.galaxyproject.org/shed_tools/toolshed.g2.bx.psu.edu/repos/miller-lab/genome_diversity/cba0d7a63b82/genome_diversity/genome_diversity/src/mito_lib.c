// mito_data.c -- shared procedures to read SNP and coverage file for
// mitogenomes

#include "lib.h"
#include "mito_lib.h"

#define ADHOC

// get the adequately covered intervals for each specified individual;
// merge adjacent intervals
void get_intervals(char *filename) {
	FILE *fp = ckopen(filename, "r");
	char buf[500], name[100];
	struct interv *p, *new;
	int i, b, e, cov;

	while (fgets(buf, 500, fp)) {
		if (sscanf(buf, "%*s %d %d %d %s", &b, &e, &cov, name) != 4)
			fatalf("interval: %s", buf);
		if (cov < min_cov)
			continue;
		for (i = 0; i < nM && !same_string(M[i].name, name); ++i)
			;
		if (i == nM)
			continue;
		if (M[i].last != NULL && M[i].last->e == b) {
			// merge with adjacent interval
			M[i].last->e = e;
			continue;
		}
		new = ckalloc(sizeof(*new));
		new->b = b;
		new->e = e;
		new->next = NULL;
		if ((p = M[i].last) == NULL)
			M[i].intervals = new;
		else 
			p->next = new;
		M[i].last = new;
	}
	fclose(fp);
/*
	for (i = 0; i < nM; ++i) {
		printf("%s:", M[i].name);
		for (p = M[i].intervals; p; p = p->next)
			printf(" %d-%d", p->b, p->e);
		putchar('\n');
	}
*/
}

// get the SNPs; for each SNP set the array of (first characters from)
// genotypes of the specified samples (individuals)
int get_variants(char *filename, struct snp *S, int refcol) {
	FILE *fp = ckopen(filename, "r");
	char buf[5000], *s, *f[101], *z = " \t\n\0";
	int i, n, c;

	for (i = 0; i <= 100; ++i)
		f[i] = NULL;
	for (n = 0; fgets(buf, 500, fp); ++n) {
		if (buf[0] == '#') {
			--n;
			continue;
		}
		if (n >= MAX_SNP)
			fatal("too many SNPs");
		if (sscanf(buf, "%*s %d", &(S[n].pos)) != 1)
			fatalf("pos : %s", buf);
		S[n].g = ckalloc((nM+1)*(sizeof(char)));
		S[n].g[nM] = '\0';
		for (i = 0; i <= 100; ++i)
			if (f[i] != NULL)
				free(f[i]);
		for (i = 1, s = strtok(buf, z); s; s = strtok(NULL, z), ++i)
			f[i] = copy_string(s);
		for (i = 0; i < nM; ++i) {
			// genotype is 2 columns past the individual's 1st
			// column
			// AD HOC RULE: IF THERE IS ONE READ OF EACH ALLELE,
			// THEN IGNORE THE SNP.
			c = M[i].col;
			if (refcol == 0)
				S[n].g[i] = f[c+2][0];
			else if (same_string(f[refcol+2], f[c+2]))
				S[n].g[i] = '2';
			else
				S[n].g[i] = '0';
#ifdef ADHOC
			if (same_string(f[c], "1") &&
			    same_string(f[c+1], "1"))
				S[n].g[i] = '-';
#endif
		}
	}
	fclose(fp);
	return n;
}
