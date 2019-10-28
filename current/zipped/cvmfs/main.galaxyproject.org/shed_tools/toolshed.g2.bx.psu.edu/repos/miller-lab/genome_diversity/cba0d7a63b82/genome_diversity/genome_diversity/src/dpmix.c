/* dpmix -- admixture using dynamic programming
*
*    argv{1] = a Galaxy SNP table. For each of several individuals, the table
*              has four columns (#A, #B, genotype, quality) -- SNPs on the same
*	       chromosome must appear together, and in order of position
*    argv[2] = column with the chromosome name (position is the next column)
*    argv[3] = "all" or e.g., "chr20"
*    argv[4] = 1 if ancestral allele frequencies are estimated from SAMtools
*		genotypes; 0 means use read-coverage data.
*    argv[5] = 1 to add logarithms of probabilities, allowing unobserve alleles,
*	       0 to simply add probabilities
*    argv[6] = switch penalty (>= 0)
*    argv[7] = file giving heterochromatic intervals ('-' means that no file is
*	       given)
*    argv[8] = file name for additional output
*    argv[9], argv[10], ...,  have the form "13:1:Peter", "13:2:Paul" or
*	      "13:0:Mary", meaning that the 13th and 14th columns (base 1)
*	      give the allele counts for an individual that is in ancestral
*	      population 1, ancestral population 2, or is a potentially admixed
*	      individual, resp.

What it does on Galaxy
The user specifies two "ancestral" populations (i.e., sources for chromosomes) and a set of potentially admixed individuals, and chooses between the sequence coverage or the estimated genotypes to measure the similarity of genomic intervals in admixed individuals to the two classes of ancestral chromosomes. The user also picks a "switch penalty", typically between 10 and 100. For each potentially admixed individual, the program divides the genome into three "genotypes": (0) homozygous for the second ancestral population (i.e., both chromosomes from that population), (1) heterozygous, or (2) homozygous for the second ancestral population. Parts of a reference chromosome that are labeled as "heterochromatic" are given the non-genotype, 3. Smaller values of the switch penalty (corresponding to more ancient admixture events) generally lead to the reconstruction of more frequent changes between genotypes.
*/

#include "lib.h"
//#include <math.h>

// maximum length of a line from the table
#define MOST 50000

// we create a linked list of "events" on a chromosome -- mostly SNPs, but
// also ends of hetorochomatic intervals
struct snp {
	double F1, F2;	// reference allele frequencies in the two populations
	int pos, *g,	// position and an array of admixed genotypes
	  type;		// 0 = SNP, 1 = start of het. interval, 2 = end
	struct snp *prev;	// we keep the list in order of decreasing pos
} *last;

// array of potentially admixed individuals
struct admixed {
	char *name;
	int gcol, ge20, gt02;
	long long x[4];		// number of reference bp in each state
} A[MOST];

// information about "ancestral" individuals, namely column and population
struct ances {
	int col, pop;
	char *name;
} C[MOST];

// heterochromatic intervals
struct het {
	char *chr;
	int b, e;
} H[MOST];

// global variables
int *B[4],	// backpointer to state at the previous SNP (or event)
    *P;		// chromosome position
int nH, nI, nG, genotypes, nsnp, debug, chr_col, logs;
char this_chr[100];
double switch_penalty;
char buf[MOST], *status;
FILE *fp, *out;

// probability of producing genotype g in admixture state s
// given reference allele frequencies f1 and f2 in the ancestral populations
double score (double f1, double f2, int g, int s) {
	double p;

	if (s == 2) { // homozygous for the first ancestral population
		if (g == 2)
			p = f1*f1;
		else if (g == 0)
			p = (1.0-f1)*(1.0-f1);
		else
			p = 2.0*f1*(1.0-f1);
	} else if (s == 0) { // homozygous for the second ancestral population
		if (g == 2)
			p = f2*f2;
		else if (g == 0)
			p = (1.0-f2)*(1.0-f2);
		else
			p = 2.0*f2*(1.0-f2);
	} else { // one chromosome from each ancestral population
		if (s != 1)
			fatalf("bad state %d", s);
		if (g == 2)
			p = f1*f2;
		else if (g == 0)
			p = (1.0-f1)*(1.0-f2);
		else
			p = f1*(1.0-f2)  + (1.0-f1)*f2;
	}
	
	if (p < 0.0)
		fatalf("%f %f %d %d => %f", f1, f2, g, s, p);
	if (!logs)
		return p;
#ifdef NEVER
	if (p == 0.0)
		return -5.0;
	p = log(p);
	if (p < -5.0)
		p = -5.0;
	return p;
#endif
	fatal("dpmix: cannot happen");
}

char *get_chr_name() {
	static char tmp[MOST];
	char *s, *z = "\t\n";
	int i = chr_col;

	strcpy(tmp, buf);
	s = strtok(tmp, z);
	while (--i > 0)
		s = strtok(NULL, z);
	return s;
}

/* Process the a-th potentially admixed individual.
*  We think of a graph with nodes (event, state) for each event (SNP or
*  end-point of a heterochromatic interval on the current chromosome) and state
*  = 0, 1, 2, 3 (corresponding to genotypes 0, 1, and 2, plus 3 =
*  heterochromatin); for events other than the last one, there are edges from
*  each (event, state) to (event+1, k) for 0 <= k <= 3. An edge (event, j) to
*  (event+1, k) has penalty 0 if j = k and penalty switch_penalty otherwise.
*  The bonus at SNP node (event, state) for 0 <= state <= 2 is the probability
*  of generating the genotype observed in the a-th potentially admixed
*  individual given the allele frequences in the two ancestral populations and
*  the assumed admixture state in this region of the chromosome. The score of a
*  path is the sum of the node bonuses minus the sum of the edge penalties.
*
*  Working backwards through the events, we compute the maximum path score,
*  from[state], from (event,state) back to the closest admixed interval.
*  To force paths to reach state 3 at an event signalling the start of a
*  heterochromatic interval (type = 1), but to avoid state 3 at other events,
*  we assign huge but arbitrary negative scores (see "avoid", below).
*  At (event,state), B[event][state] is the backpointer to the state at
*  event+1 on an optimal path. Finally, we follow backpointers to partition
*  the chromosome into admixture states.
*/
void one_admix(int a) {
	int i, j, m, state, prev_pos, b;
	double from[4], f[4], ff[4], avoid = -1000000.0;
	struct snp *p;

	// from[i] = highest score of a path from the current event
	// (usually a SNP) to the next (to the right) heterochromatic interval
	// or the end of the chromosome. The score of the path is the sum of
	// SNP scores minus (switch_penalty times number of state switches). 
	// We assume that the last two event on the chromosome are the start
	// and end of a heterochromatic interval (possibly of length 0)/
	for (i = 0; i < 4; ++i)
		from[i] = 0;
	for (i = nsnp-1, p = last; i >= 0 && p != NULL; --i, p = p->prev) {
		if (p->type == 0 && p->g[a] == -1) { // no genotype
			// no state change at this SNP
			for (state = 0; state < 4; ++state)
				B[state][i] = state;
			continue;
		}
		
		for (state = 0; state < 4; ++state) {
			// find highest path-score from this event onward
			for (m = j = 0; j < 4; ++j) {
				f[j] = from[j];
				if (j != state)
					f[j] -= switch_penalty;
				//if (abs(j-state) == 2)
					//from[j] -= switch_penalty;
				if (f[j] > f[m])
					m = j;
			}
			B[state][i] = m;
			ff[state] = f[m];
			if (state < 3 && p->type == 0)
				ff[state] +=
				    score(p->F1, p->F2, p->g[a], state);
		}
		if (p->type == 1) {
			// start of heterochomatic interval. Force paths
			// reaching this point to go through state 3
			from[3] = 0;
			from[0] = from[1] = from[2] = avoid;
		} else {
			for (j = 0; j < 3; ++j)
				from[j] = ff[j];
			from[3] = avoid;
		}
		if (debug)
			fprintf(stderr, "%d: %f(%d) %f(%d) %f(%d) %f(%d)\n",
			  i, from[0], B[0][i], from[1], B[1][i], from[2],
			  B[2][i], from[3], B[3][i]);
	}

	// find the best initial state
	for (state = 0, j = 1; j < 4; ++j)
		if (from[j] > from[state])
			state = j;

	// trace back to find the switch points
	// A[a].x[state] records the total length of intervals in each state
	for (prev_pos = i = 0; i < nsnp; ++i) {
		if ((b = B[state][i]) != state) {
			if (prev_pos < P[i+1]-1)
				printf("%s\t%d\t%d\t%d\t%s\n",
				  this_chr, prev_pos, P[i+1], state, A[a].name);
			A[a].x[state] += (P[i+1]-prev_pos);
			prev_pos = P[i+1];
			state = b;
		}
	}
} 

// Add a heterochromatic interval to the SNP list, where type = 1 signifies
// the start of the interval, 2 signifies the end.
void add_het(int b, int type) {
	struct snp *new = ckalloc(sizeof(struct snp));
	int i;

	new->F1 = new->F2 = 0.0;
	new->pos = b;
	new->type = type;
	new->g = ckalloc(nG*sizeof(int));
	for (i = 0; i < nG; ++i)
		new->g[i] = 0;
	new->prev = last;
	last = new;
}

/* Process one chromosome. Read the SNPs on the chromosome (the first one is
*  already in the buf). Boil each SNP down to the contents of a SNP entry
*  (pos, F1, F2, g[]) and put it in the linked list. Also, intersperse the
*  "events" corresponding to the start and end of a heterochromatic interval.
*  Then call the dynamic-programming routine for each potentially admixed
*  individual.
*/
void one_chr() {
	char *s, *z = "\t\n";
	int X[MOST], n, i, g, A1, B1, A2, B2, a, do_read, p, pos, het;
	struct snp *new;
	double F1, F2;

	strcpy(this_chr, get_chr_name());
	nsnp = 0;
	last = NULL;
	// advance to this chromosome in the list of heterochromatic intervals
	for (het = 0; het < nH && !same_string(this_chr, H[het].chr); ++het)
		;
	// loop over the SNPs on the current chromosome
	for (do_read = 0; ; do_read = 1) {
		if (do_read && (status = fgets(buf, MOST, fp)) == NULL)
			break; 
		if (!same_string(get_chr_name(), this_chr))
			break;
		
		// set X[i] = atoi(i-th word of buf), i is base 1
		for (i = 1, s = strtok(buf, z); s != NULL;
		  ++i, s = strtok(NULL, z))
			X[i] = atoi(s);

		// insert events (pseudo-SNPs) for heterochomatin intervals
		// coming before the SNP
		pos = X[chr_col+1];
		while (het < nH && same_string(this_chr, H[het].chr) &&
		   H[het].b < pos) {
			add_het(H[het].b, 1);
			add_het(H[het].e, 2);
			nsnp+= 2;
			++het;
		}
			
		// should we discard this SNP?
		if (pos == -1)	// SNP not mapped to the reference
			continue;
/*
		for (i = 0; i < nG && X[A[i].gcol] >= 0; ++i)
			;
		if (i < nG)	// genotype of admixed individual not called
			continue;
*/

		// add SNP to a "backward pointing" linked list, recording the
		// major allele frequencies in the two reference populations
		// and genotypes in the potential admixed individuals
		for (i = A1 = B1 = A2 = B2 = 0; i < nI; ++i) {
			n = C[i].col;
			p = C[i].pop;
			if (genotypes) {
				g = X[n+2];
				if (g == -1)
					continue;
				if (g < 0 || g > 2)
					fatalf("invalid genotype %d", g);
				if (p == 1) {
					A1 += g;
					B1 += (2 - g);
				} else if (p == 2) {
					A2 += g;
					B2 += (2 - g);
				}
			} else {	// use read counts
				if (p == 1) {
					A1 += X[n];
					B1 += X[n+1];
				} else if (p == 2) {
					A2 += X[n];
					B2 += X[n+1];
				}
			}
		}
		if (A1+B1 == 0 || A2+B2 == 0)
			continue;
		++nsnp;
		new = ckalloc(sizeof(struct snp));
		new->pos = X[chr_col+1];
		new->F1 = F1 = (double)A1/(double)(A1+B1);
		new->F2 = F2 = (double)A2/(double)(A2+B2);
		new->type = 0;
		new->g = ckalloc(nG*sizeof(int));
		for (i = 0; i < nG; ++i) {
			g = new->g[i] = X[A[i].gcol];
			if (score(F1, F2, g, 2) >= score(F1, F2, g, 0))
				A[i].ge20++;
			else 
				A[i].gt02++;
		}
		if (F1 < 0.0 || F1 > 1.0)
			fatalf("F1 = %f (A1 = %d, B1 = %d) at snp %d",
			  F1, A1, B1, nsnp);
		if (F2 < 0.0 || F2 > 1.0)
			fatalf("F2 = %f (A2 = %d, B2 = %d) at snp %d",
			  F2, A2, B2, nsnp);
		new->prev = last;
		last = new;
	}
	// insert heterochomatin intervals that follow all SN
	while (het < nH && same_string(this_chr, H[het].chr)) {
		add_het(H[het].b, 1);
		add_het(H[het].e, 2);
		nsnp += 2;
		++het;
	}
/*
printf("nsnp = %d\n", nsnp);
for (i = nsnp-1, new = last; i >= 0 && new != NULL; --i, new = new->prev) {
printf("%d %d ", new->pos, new->type);
printf("%g %g ", new->F1, new->F2);
for (a = 0; a < nG; ++a)
printf("%d", new->g[a]);
putchar('\n');
}
//exit(0);
printf("\nbacktrace\n");
*/

	// allocate arrays for the DP analysis
	P = ckalloc(nsnp*sizeof(int));	// position of each event
	for (i = nsnp-1, new = last; i >= 0 && new != NULL;
	     --i, new = new->prev)
		P[i] = new->pos;

	for (i = 0; i < 4; ++i) {	// space for back-pointers
		B[i] = ckalloc((nsnp+1)*sizeof(int));
		B[i][nsnp] = 0;
	}
	
	// loop over possibly admixed individuals
	for (a = 0; a < nG; ++a)
		one_admix(a);

	// free the allocated storage
	while (last != NULL) {
		new = last;
		last = last->prev;
		free(new->g);
		free(new);
	}
	free(P);
	for (i = 0; i < 4; ++i)
		free(B[i]);
}

int main(int argc, char **argv) {
	int n, i, j, k, saw[3];
	long long het_len, ref_len;
	float N;
	char nam[100], *chr;

	if (argc < 9)
		fatal("args: table chr-col chr data-source logs switch heterochrom outfile n:1:name1 m:2:name2 ...");
	if (same_string(argv[argc-1], "debug")) {
		debug = 1;
		--argc;
	}

	// handle command-line arguments
	chr_col = atoi(argv[2]);
	chr = argv[3];
	genotypes = atoi(argv[4]);

	logs = atoi(argv[5]);
	if (logs)
		fatal("logarithms of probabilities -- under development");
	//if (logs) switch_penalty = log(switch_penalty);

	switch_penalty = atof(argv[6]);
	if (switch_penalty < 0.0)
		fatal("negative switch penalty");
	out = ckopen(argv[8], "w");

	het_len = ref_len = 0;
	if (!same_string(argv[7], "-")) {
		fp = ckopen(argv[7], "r");
		while (fgets(buf, MOST, fp)) {
			if (nH >= MOST)
				fatal("Too many heterochromatic intervals");
			if (sscanf(buf, "%s %d %d", nam, &i, &j) != 3)
				fatalf("gagging: %s", buf);
			H[nH].chr = copy_string(nam);
			H[nH].b = i;
			H[nH].e = j;
			// assumes last event per chrom. is a het. interval
			if (nH > 0 && !same_string(nam, H[nH-1].chr))
				ref_len += j;
			het_len += (j - i);
			++nH;
		}
		fclose(fp);
	}
	ref_len += H[nH-1].e;

	// populations must be disjoint
	saw[1] = saw[2] = 0;
	for (i = 9; i < argc; ++i) {
		if (sscanf(argv[i], "%d:%d:%s", &j, &k, nam) != 3)
			fatalf("not like 13:2:fred : %s", argv[i]);
		if (k < 0 || k > 2)
			fatalf("not population 0, 1 or 2: %s", argv[i]);
		saw[k] = 1;

		// seen this individual (i.e., column) before??
		for (n = 0; n < nI && C[n].col != j; ++n)
			;
		if (n < nI)
			fatal("populations are not disjoint");
		if (k == 0) {	// admixed individual
			if (nG >= MOST)
				fatal("Too many admixed individuals");
			A[nG].name = copy_string(nam);
			A[nG++].gcol = j+2;
		} else {	// in an ancestral population
			if (nI >= MOST)
				fatal("Too many ancestral individuals");
			C[nI].col = j;
			C[nI].pop = k;
			C[nI++].name = copy_string(nam);
		}
	}
	if (saw[0] == 0)
		fatal("no admixed individual is specified");
	if (saw[1] == 0)
		fatal("first reference population is empty");
	if (saw[2] == 0)
		fatal("second reference population is empty");

	// start the output file of text
	for (k = 1; k <= 2; ++k) {
		fprintf(out, "state %d agrees with:", k == 1 ? 2 : 0);
		for (i = 0; i < nI; ++i)
			if (C[i].pop == k)
				fprintf(out, " %s", C[i].name);
		putc('\n', out);
	}
	putc('\n', out);

	fp = ckopen(argv[1], "r");
	while ((status = fgets(buf, MOST, fp)) != NULL && buf[0] == '#')
		;
	if (same_string(chr, "all"))
		while (status != NULL)
			one_chr();
	else {	// skip to the specified chromosome
		while (!same_string(chr, get_chr_name()) &&
		       (status = fgets(buf, MOST, fp)) != NULL)
			;
		if (status != NULL)
			one_chr();
	}
	for (i = 0; i < nG; ++i) {
		fprintf(out,
		  "%s: %d SNPs where state 2 is at least as likely as state 0\n",
		  A[i].name, A[i].ge20);
		fprintf(out,
		  "%s: %d SNPs where state 0 is more likely than state 2\n\n",
		  A[i].name, A[i].gt02);
	}
	// write fractions in each state to the output text file

	if (ref_len)
		fprintf(out,
		  "%lld of %lld reference bp (%1.1f%%) are heterochromatin\n\n",
		  het_len, ref_len, 100.0*(float)het_len/(float)ref_len);

	for (i = 0; i < nG; ++i) {
		N = (float)(A[i].x[0] + A[i].x[1] + A[i].x[2])/100.0;
		fprintf(out, "%s: 0 = %1.1f%%, 1 = %1.1f%%, 2 = %1.1f%%\n",
		  A[i].name, (float)A[i].x[0]/N, (float)A[i].x[1]/N,
		  (float)A[i].x[2]/N); 
	}

	return 0;
}
