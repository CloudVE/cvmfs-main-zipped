/* mk_Ji -- prepare data for drawing a picture of mitogenome data
*
*  argv[1] -- SNP table for the mitogenome
*
*  argv[2] -- indel table for the mitogenome
*
*  argv[3] -- coverage table: file of intervals with lines like:

P.ingens-mt     175     194     6       9-M-352

*  giving genome name, start postion (base-0), end position (half open),
*  coverage and sample name.
*
*  argv[4] -- annotation file like

P.ingens-mt     0       70      tRNA    +       tRNA-Phe
P.ingens-mt     70      1030    rRNA    +       12S
P.ingens-mt     1030    1100    tRNA    +       tRNA-Val
P.ingens-mt     1101    2680    CDS     +       rRNA
P.ingens-mt     1101    2680    rRNA    +       16S
P.ingens-mt     2680    2755    tRNA    +       tRNA-Leu
P.ingens-mt     2758    3713    CDS     +       ND1
...
P.ingens-mt     15484   16910   D-loop  +       D-loop

*  argv[5] -- the minimum coverage. Intervals of lower coverage are ignored.
*
*  argv[6] -- either the string "default" or the name of an individual
*
*  argv[7], argv[8], ... column:name pairs like "9:sam".
*  
*  Also, if the last argument is "debug", then much output sent to stderr.
*/

#include "lib.h"
#include "mito_lib.h"

int ref_len;

// read gene annotation, change "CDS" to "gene", print for the drawing tool,
// print lines showing the genome name and length (last annotated position).
void get_genes(char *filename) {
	FILE *fp = ckopen(filename, "r");
	char buf[500], ref[100], type[100], name[100], *t;
	int b, e;

	while (fgets(buf, 500, fp)) {
		if (sscanf(buf, "%s %d %d %s %*s %s",
		    ref, &b, &e, type, name) != 5)
			fatalf("gag: %s", buf);
		t = (same_string(type, "CDS") ? "gene" : type);
		// print the Genome Annotation line
		printf("@GA=%d:%d:%s:%s\n", b, e, name, t);
	}
	printf("@GL=%d\n", ref_len = e);
	printf("@GN=%s\n", ref);
}

// print items that are adequately covered
void visible(int i, struct snp *S, int nS, char *s) {
	struct interv *t;
	int j;

	for (j = 0, t = M[i].intervals; j < nS; ++j) {
		while (t && t->e <= S[j].pos)
			t = t->next;
		if (t && t->b <= S[j].pos && S[j].g[i] == '0')
			printf(" %s%d", s, S[j].pos);
	}
}

int main(int argc, char **argv) {
	struct interv *t;
	int i, nS, nI, last_e, refcol;
	char *a, *s;

	if (argc > 7 && same_string(argv[argc-1], "debug")) {
		--argc;
		debug = 1;
	}

	if (argc < 7)
		fatal("args: snps indels intervals genes min_cov ref 9:sam 13:judy ... ");

	// store sample names and start positions (argv[6], argv[7], ...)
	for (nM = 0, i = 7; i < argc; ++nM, ++i) {
		if (nM >= MAX_SAMPLE)
			fatalf("Too many mitogenomes");
		if ((s = strchr(a = argv[i], ':')) == NULL)
			fatalf("colon: %s", a);
		M[nM].col = atoi(a);
		M[nM].name = copy_string(s+1);
	}
	if (same_string(argv[6], "default"))
		refcol = 0;
	else {
		for (i = 0; i < nM && !same_string(argv[6], M[i].name); ++i)
			;
		if (i == nM)
			fatalf("improper reference name '%s'", argv[6]);
		refcol = M[i].col;
		// fprintf(stderr, "refcol = %d\n", refcol);
	}

	// read annotation and annotate in the file for drawing
	get_genes(argv[4]);

	// record color information
	printf("@CL=rRNA:#EF8A62\n@CL=tRNA:#31A354\n@CL=gene:#B2182B\n");
	printf("@CL=missing:#67A9CF:L\n@CL=indel:#2166AC\n@CL=special:#000000\n");

	min_cov = atoi(argv[5]);

	// store the coverage data
	get_intervals(argv[3]);

	if (debug) {
		for (i = 0; i < nM; ++i) {
			fprintf(stderr, ">%s\n", M[i].name);
			for (t = M[i].intervals; t; t = t->next)
				fprintf(stderr, "%d\t%d\n", t->b, t->e);
		}
		putc('\n', stderr);
	}

	// record the variants
	nS = get_variants(argv[1], S, refcol);
	nI = get_variants(argv[2], I, refcol);

	// report the information for each sample
	for (i = 0; i < nM; ++i) {
		printf("%s", M[i].name);
		visible(i, S, nS, "");
		visible(i, I, nI, "indel=");
		last_e = 0;
		for (t = M[i].intervals; t; t = t->next) {
			if (last_e < t->b)
				printf(" missing=%d:%d", last_e, t->b);
			last_e = t->e;
		}
		if (last_e < ref_len)
			printf(" missing=%d:%d", last_e, ref_len);
		putchar('\n');
	}

	return 0;
}
