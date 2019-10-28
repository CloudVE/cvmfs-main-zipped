/* mt_pi -- compute the diversity measure pi for mitochondrial genomes
*  [SHOULD I OPTIONALLY INCLUDE INDELS?]
*
*  argv[1] -- SNP table for the mitogenome
*
*  argv[2] -- file of intervals with lines like:

P.ingens-mt     175     194     6       9-M-352

*  giving genome name, start postion (base-0), end position (half open),
*  coverage and sample name.
*
*  argv[3] -- the minimum coverage. Intervals of lower coverage are ignored.
*
*  argv[4], argv[5], ... column:name pairs like "9:sam".
*  
*  Also, if the last argument is "debug", then much output sent to stderr, if it
*  is "debug2", then the coverage and difference-count for each mitogenome-pair
*  is sent to stderr.
*/

#include "lib.h"
#include "mito_lib.h"

int debug2;

// for a pair of samples, determine how much of the reference is in the
// adequately covered segments for both, and count the number of SNPs in those
// shared regions where they differ
// PUTATIVE HETEROPLASMIES ARE IGNORED
float pair(int i, int j, int nS) {
	int covered, B, E, diffs, k;
	struct interv *p = M[i].intervals, *q = M[j].intervals;
	char x, y;

	// k scans the SNPs
	covered = diffs = k = 0;
	while (p && q) {
		if (debug)
			fprintf(stderr, "trying %d-%d and %d-%d\n",
			  p->b, p->e, q->b, q->e);
		// take the intersection of the two well-covered intervals
		B = MAX(p->b, q->b);
		E = MIN(p->e, q->e);
		if (B < E) {
			if (debug)
				fprintf(stderr, "  covered %d\n", E-B);
			covered += (E - B);
			for ( ; k < nS && S[k].pos < E; ++k) {
				if (S[k].pos >= B) {
					x = S[k].g[i];
					y = S[k].g[j];
					if (debug)
						fprintf(stderr,
						  "   SNP %c %c at %d\n",
						  x, y, S[k].pos);
/*
					if (x == '-' || y == '-')
						fatalf("SNP at %d missing geno",
						  S[k].pos);
*/
/*
					if (x == '1' || y == '1')
						continue;
*/
					if (x != y) {
						++diffs;
						if (debug)
							fprintf(stderr,
							  "\tdiff at %d\n",
							  S[k].pos);
					}
				}
			}
		}
		// go to next adequately covered interval(s)
		if (p->e < q->e)
			p = p->next;
		else if (p->e > q->e)
			q = q->next;
		else {
			p = p->next;
			q = q->next;
		}
	}
	if (debug2)
		fprintf(stderr, "cov(%s,%s) = %d, diffs = %d\n",
		  M[i].name, M[j].name, covered, diffs);
/*
	if (covered == 0)
		fatalf("coverage threshold is too high for %s and %s",
		  M/[i].name, M[j].name);
*/
	if (covered == 0)
		return -1.0;
	return (float)diffs/(float)covered;
}

int main(int argc, char **argv) {
	struct interv *t;
	int i, j, nS, good_pairs, bad_pairs;
	char *a, *s;
	float tot, pr;

	if (argc > 4 && same_string(argv[argc-1], "debug")) {
		--argc;
		debug = debug2 = 1;
	} else if (argc > 4 && same_string(argv[argc-1], "debug2")) {
		--argc;
		debug2 = 1;
	}

	if (argc < 5)
		fatal("args: snps intervals min_cov 9:sam 13:judy ... ");
	// store sample names and start positions (argv[4], argv[5], ...)
	for (nM = 0, i = 4; i < argc; ++nM, ++i) {
		if (nM >= MAX_SAMPLE)
			fatalf("Too many mitogenomes");
		if ((s = strchr(a = argv[i], ':')) == NULL)
			fatalf("colon: %s", a);
		M[nM].col = atoi(a);
		M[nM].name = copy_string(s+1);
	}
	min_cov = atoi(argv[3]);
	get_intervals(argv[2]);

	if (debug) {
		for (i = 0; i < nM; ++i) {
			fprintf(stderr, ">%s\n", M[i].name);
			for (t = M[i].intervals; t; t = t->next)
				fprintf(stderr, "%d\t%d\n", t->b, t->e);
		}
		putc('\n', stderr);
	}

	// record the SNPs
	nS = get_variants(argv[1], S, 0);

	if (debug) {
		for (i = 0; i < nS; ++i)
			fprintf(stderr, "%d %s\n", S[i].pos, S[i].g);
		putc('\n', stderr);
	}

	// record the total rate of diversity, over all pairs of individuals
	// having overlapping well-covered intervals
	good_pairs = bad_pairs = 0;
	for (i = 0, tot = 0.0; i < nM; ++i) {
		for (j = i+1; j < nM; ++j) {
			pr = pair(i, j, nS);
			if (pr >= 0.0) {
				++good_pairs;
				tot += pr;
			} else
				++bad_pairs;
		}
	}
	printf("pi = %5.5f\n", tot/(float)good_pairs);
	if (bad_pairs > 0)
		printf("%d of %d pairs had no sequenced regions in common\n",
		  bad_pairs, bad_pairs + good_pairs);

	return 0;
}
