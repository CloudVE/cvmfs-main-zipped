/* dpmix -- admixture using dynamic programming; 2 or 3 source populations
*
*    argv{1] = a Galaxy SNP table. For each individual, the table may have
*              four columns (#A, #B, genotype, quality), or only one
*	       (genotype). SNPs on the same chromosome must appear together,
*	       and in order of position
*    argv[2] = column with the chromosome name (position is the next column)
*    argv[3] = "all" or e.g., "chr20"
*    argv[4] = 1 if source-pop allele frequencies are estimated from genotypes;
*	       0 means use read-coverage data.
*    argv[5] = 1 to add logarithms of probabilities;
*	       0 to simply add probabilities
*    argv[6] = switch penalty (>= 0)
*    argv[7] = file giving heterochromatic intervals ('-' = no file is given)
*    argv[8] = file name for additional output
*    argv[9], argv[10], ... =  like "13:1:Peter", "13:2:Paul", "13:3:Sam"
*	       or "13:0:Mary", meaning that the 13th and 14th columns (base 1)
*	       give the allele counts for an individual that is in source
*	       population 1, source population 2, source population 3,
*	       or is a potentially admixed individual, resp. Even when only the
*	       genotype is given, it is found in column 15.
*    optional last arguments = "debug"

What it does on Galaxy
The user specifies two or three source populations (i.e., sources for chromosomes) and a set of potentially admixed individuals, and chooses between the sequence coverage or the estimated genotypes to measure the similarity of genomic intervals in admixed individuals to the source chromosomes. The user also specifies a "switch penalty", controlling the strength of evidence needed to switch between source populations as the the program scans along a chromosome.

Choice of an appropriate switch penalty depends on the number of SNPs and the time since the admixture event(s), and it is not always easy to determine a good setting ahead of time. One approach is to try, say, 10.0, and look at the resulting picture. If the admixture blocks are shorter than anticipated, based on expectations about the number of generations since the admixture event(s) (i.e., the implied admixture events are too ancient), then try a larger value. Conversely, if the blocks are longer than anticipated (i.e., the implied timing of admixture is too recent), then lower the penalty. In any case, it is prudent to base conclusions on results consistently obtained with a range of switch penalities.

If there are 3 source populatons, then for each potentially admixed individual the program divides each potentially admixed genome into six "genotypes":

(1) homozygous for the first source population (i.e., both chromosomes from that population),
(2) homozygous for the second source population,
(3) homozygous for the third source population,
(4) heterozygous for the first and second populations (i.e., one chromosome from each),
(5) heterozygous for the first and third populations, or
(6) heterozygous for the second and third populations.

Parts of a reference chromosome that are labeled as "heterochromatic" are given the "non-genotype" 0. With two source populations, there are only three possible "genotypes".

There are two Galaxy history-items generated. The first is a tabular dataset with chromosome, start, stop, and pairs of columns containing the "genotypes" as described above and the label from the admixed individual. It is used to draw a pciture and can be consulted for details. The second item is a composite of (1) a link to a pdf which graphically shows the inferred source population for each potentially admixed individual along each chromosome, and (2) a link to a text file that describes the run and summarizes the extent of predicted admixture in each individual.

For certain genome assemblies, Galaxy stores a "heterochromatin file", which specifies the heterochromatic intervals on each chromosome, as well as telling the chromosome lengths. For other genome assemblies, the user can supply this information. Here are the first few lines of a heterochromatin file for human assembly hg19:
chr1 121485434 142535434
chr1 249250621 249250621
chr2 90545103 95326171
chr2 243199373 243199373
This gives the centromeres and lengths for chromosomes 1 and 2. Chromosomal addresses begin at 0, and the stated end position is 1 beyond the last address. For instance the interval described by the second line has length 0; it tells Galaxy that chromosome 1 is 249,250,621 bp in length. The entries for an acrocentric chromosome looks like the following:
chr22 0 16050000
chr22 51304566 51304566
The file must be ordered by chromosome name (starting with autosomes), and by position within a chromosome.
*/

#include "lib.h"
#include <math.h>

// maximum length of a line from the table
#define MOST 50000

// we create a linked list of "events" on a chromosome -- mostly SNPs, but
// also ends of hetorochomatic intervals
struct snp {
	double F1, F2, F3; // ref. allele frequencies in the three populations
	int pos, *g,	// position; genotypes of admixed individuals
	  type;		// 0 = SNP, 1 = start of het. interval, 2 = end
	struct snp *prev;	// we keep the list in order of decreasing pos
} *last;

// array of potentially admixed individuals
struct admixed {
	char *name;
	int gcol;
	double x[7];		// number of reference bp in each state
} A[MOST];

// information about source-population individuals
struct ances {
	int col, pop;
	char *name;
} C[MOST];

// heterochromatic intervals
struct het {
	char *chr;
	int b, e;
} H[MOST];

#define MAX_CHR_NAME 1000000

// global variables
int *B[7],	// backpointer to state at the previous SNP (or event)
    *P;		// chromosome position
int nH, nI, nG, genotypes, nsnp, debug, chr_col, logs, last_snp_pos, pop3,
 nchr_name;
char this_chr[100], *chr_name[MAX_CHR_NAME];
double switch_penalty;
char buf[MOST], *status;
FILE *fp, *out;

#define LOG_ZERO -4.0
// probability of producing genotype g in admixture state s = 1 to 6
// given source-population ref. allele frequencies f1, f2, and f3
double score (double f1, double f2, double f3, int g, int s) {
	double p;

	if (g < 0 || g > 2)
		fatalf("bad genotype %d", g);
	if (s == 1) { // homozygous for the first source population
		if (g == 2)
			p = f1*f1;
		else if (g == 0)
			p = (1.0-f1)*(1.0-f1);
		else
			p = 2.0*f1*(1.0-f1);
	} else if (s == 2) { // homozygous for the second source population
		if (g == 2)
			p = f2*f2;
		else if (g == 0)
			p = (1.0-f2)*(1.0-f2);
		else
			p = 2.0*f2*(1.0-f2);
	} else if (s == 3) { // homozygous for the third source population
		if (g == 2)
			p = f3*f3;
		else if (g == 0)
			p = (1.0-f3)*(1.0-f3);
		else
			p = 2.0*f3*(1.0-f3);
	} else if (s == 4) { // one chromosome each from source pops 1 and 2
		if (g == 2)
			p = f1*f2;
		else if (g == 0)
			p = (1.0-f1)*(1.0-f2);
		else
			p = f1*(1.0-f2)  + (1.0-f1)*f2;
	} else if (s == 5) { // one chromosome each from source pops 1 and 3
		if (g == 2)
			p = f1*f3;
		else if (g == 0)
			p = (1.0-f1)*(1.0-f3);
		else
			p = f1*(1.0-f3)  + (1.0-f1)*f3;
	} else if (s == 6) { // one chromosome each from source pops 2 and 3
		if (g == 2)
			p = f2*f3;
		else if (g == 0)
			p = (1.0-f2)*(1.0-f3);
		else
			p = f2*(1.0-f3)  + (1.0-f2)*f3;
	} else
		fatalf("bad state %d", s);
	
	if (p < 0.0)
		fatalf("%f %f %f %d %d => %f", f1, f2, f3, g, s, p);
	if (!logs)
		return p;
	if (p == 0.0)
		return LOG_ZERO;
	p = MAX(log(p), LOG_ZERO);
	return p;
	fatal("dpmix: cannot happen");
}

char *get_chr_name() {
	static char tmp[MOST];
	char *s, *z = "\t\n";
	int i = chr_col;
	static int autosome_warning = 0;

	strcpy(tmp, buf);
	s = strtok(tmp, z);
	while (--i > 0)
		s = strtok(NULL, z);
	if (!autosome_warning && strncmp(s, "chr", 3) == 0 &&
	    !isdigit(s[3])) {
		fprintf(out,
		  "WARNING: results assume diploid (non-sex) chromosomes\n\n");
		autosome_warning = 1;
	}
	return s;
}

/* Process the a-th potentially admixed individual.
*  We think of a graph with nodes (event, state) for each event (SNP or
*  end-point of a heterochromatic interval on the current chromosome) and state
*  = 0 through 7 (corresponding to genotypes 1 to 6, plus 0 =
*  heterochromatin); for events where state != 0, there are 7 edges from
*  each (event, state) to (event+1, k) for 0 <= k <= 6. An edge (event, j) to
*  (event+1, k) has penalty 0 if j = k and penalty switch_penalty otherwise.
*  The bonus at SNP node (event, state) for 1 <= state <= 6 is the probability
*  of generating the genotype observed in the a-th potentially admixed
*  individual given the allele frequences in the source populations and
*  the assumed admixture state in this region of the chromosome. The score of a
*  path is the sum of the node bonuses minus the sum of the edge penalties.
*
*  Working backwards through the events, we compute the maximum path score,
*  from[state], from (event,state) back to the closest heterochromatin interval.
*  To force paths to reach state 0 at an event signalling the start of a
*  heterochromatic interval (type = 1), but to avoid state 0 at other events,
*  we assign huge but arbitrary negative scores (see "avoid", below).
*  At (event,state), B[event][state] is the backpointer to the state at
*  event+1 on an optimal path. Finally, we follow backpointers to partition
*  the chromosome into admixture states.
*/
void one_admix(int a) {
	int i, j, m, state, prev_pos, b;
	double from[7], f[7], ff[7], avoid = -1000000.0;
	struct snp *p;

	// from[i] = highest score of a path from the current event
	// (usually a SNP) to the next (to the right) heterochromatic interval
	// or the end of the chromosome. The score of the path is the sum of
	// SNP scores minus (switch_penalty times number of state switches). 
	// We assume that the last two events on the chromosome are the start
	// and end of a heterochromatic interval (possibly of length 0)
	for (i = 0; i < 7; ++i)
		from[i] = 0;
	for (i = nsnp-1, p = last; i >= 0 && p != NULL; --i, p = p->prev) {
		if (p->type == 0 && p->g[a] == -1) { // no genotype
			// no state change at this SNP
			for (state = 0; state < 7; ++state)
				B[state][i] = state;
			continue;
		}
		
		for (state = 0; state < 7; ++state) {
			// find highest path-score from this event onward
			for (m = j = 0; j < 7; ++j) {
				f[j] = from[j];
				if (j != state)
					f[j] -= switch_penalty;
				if (f[j] > f[m])
					m = j;
			}
			B[state][i] = m;
			ff[state] = f[m];
			if (state > 0 && p->type == 0)
				ff[state] +=
				    score(p->F1, p->F2, p->F3, p->g[a], state);
		}
		if (p->type == 1) {
			// start of heterochomatic interval. Force paths
			// reaching this point to go through state 0
			from[0] = 0;
			for (j = 1; j < 7; ++j)
				from[j] = avoid;
		} else {
			for (j = 1; j < 7; ++j)
				from[j] = ff[j];
			from[0] = avoid;
		}
		if (debug) {
			fprintf(stderr, "%d:", i);
			for (j = 0; j < 7; ++j) {
				if (pop3 || j == 3 || j == 5 || j == 6)
				  fprintf(stderr, " %f(%d)", from[j], B[j][i]);
			}
			putc('\n', stderr);
		}
	}

	// find the best initial state
	for (state = 0, j = 1; j < 7; ++j)
		if (from[j] > from[state])
			state = j;

	// trace back to find the switch points
	// A[a].x[state] records the total length of intervals in each state
	for (prev_pos = i = 0; i < nsnp; ++i) {
		if ((b = B[state][i]) != state) {
			if (prev_pos < P[i+1]-1)
				printf("%s\t%d\t%d\t%d\t%s\n",
				  this_chr, prev_pos, P[i+1],
				  (state==4 && !pop3 ? 3 : state), A[a].name);
			A[a].x[state] += (double)(P[i+1]-prev_pos);
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

	new->F1 = new->F2 = new->F3 = 0.0;
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
	int X[MOST], n, i, g, A1, B1, A2, B2, A3, B3, a, do_read, p, pos, het,
	  old_pos;
	struct snp *new;
	double F1, F2, F3;

	strcpy(this_chr, get_chr_name());
	if (nchr_name == 0)
		chr_name[nchr_name++] = copy_string(this_chr);
	old_pos = nsnp = 0;
	last = NULL;
	// advance to this chromosome in the list of heterochromatic intervals
	for (het = 0; het < nH && !same_string(this_chr, H[het].chr); ++het)
		;
	// loop over the SNPs on the current chromosome
	for (do_read = 0; ; do_read = 1) {
		if (do_read && (status = fgets(buf, MOST, fp)) == NULL)
			break; 
		if (!same_string(s = get_chr_name(), this_chr)) {
			if (nchr_name >= MAX_CHR_NAME)
				fatal("Too many chromosome names");
			for (i = 0;
			     i < nchr_name && !same_string(s, chr_name[i]); ++i)
				;
			if (i < nchr_name)
				fatalf("SNVs on %s aren't together", s);
			chr_name[nchr_name++] = copy_string(s);
			break;
		}
		
		// set X[i] = atoi(i-th word of buf), i is base 1
		for (i = 1, s = strtok(buf, z); s != NULL;
		  ++i, s = strtok(NULL, z))
			X[i] = atoi(s);

		// insert events (pseudo-SNPs) for heterochomatin intervals
		// coming before the SNP
		pos = X[chr_col+1];
		if (pos <= old_pos)
			fatalf("SNP at %s %d is out of order", this_chr, pos);
		old_pos = pos;
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

		// add SNP to a "backward pointing" linked list, recording the
		// major allele frequencies in the source populations and
		// genotypes in the potential admixed individuals
		for (i = A1 = B1 = A2 = B2 = A3 = B3 = 0; i < nI; ++i) {
			n = C[i].col;
			p = C[i].pop;
			if (genotypes) {
				g = X[n+2];
				if (g == -1)
					continue;
				if (g < 0 || g > 2)
					fatalf("invalid genotype %d in column %d, pos %d", g, n+2, X[2]);
				if (p == 1) {
					A1 += g;
					B1 += (2 - g);
				} else if (p == 2) {
					A2 += g;
					B2 += (2 - g);
				} else if (p == 3) {
					A3 += g;
					B3 += (2 - g);
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
		new->F3 = F3 = (double)A3/(double)(A3+B3);
		new->type = 0;
		new->g = ckalloc(nG*sizeof(int));
		for (i = 0; i < nG; ++i)
			g = new->g[i] = X[A[i].gcol];
		if (F1 < 0.0 || F1 > 1.0)
			fatalf("F1 = %f (A1 = %d, B1 = %d) at snp %d",
			  F1, A1, B1, nsnp);
		if (F2 < 0.0 || F2 > 1.0)
			fatalf("F2 = %f (A2 = %d, B2 = %d) at snp %d",
			  F2, A2, B2, nsnp);
		if (F3 < 0.0 || F3 > 1.0)
			fatalf("F3 = %f (A2 = %d, B2 = %d) at snp %d",
			  F3, A3, B3, nsnp);
		new->prev = last;
		last = new;
	}
		
	// insert heterochomatin intervals that follow all SNPs
	while (het < nH && same_string(this_chr, H[het].chr)) {
		add_het(H[het].b, 1);
		add_het(H[het].e, 2);
		nsnp += 2;
		++het;
	}
	// make sure the picture is drawn to at least the last SNP
	if (last->type == 0) {	
		i = last->pos + 1;
		add_het(i, 1);
		add_het(i, 2);
		nsnp += 2;
	}

	// allocate arrays for the DP analysis
	P = ckalloc(nsnp*sizeof(int));	// position of each event
	for (i = nsnp-1, new = last; i >= 0 && new != NULL;
	     --i, new = new->prev)
		P[i] = new->pos;

	for (i = 0; i < 7; ++i) {	// space for back-pointers
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
	for (i = 0; i < 7; ++i)
		free(B[i]);
}

int main(int argc, char **argv) {
	int n, i, j, k, saw[4];
	long long het_len, ref_len;
	double N, tot[7], keep[7], xx, yy;
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
			if (nH > 0 && !same_string(nam, H[nH-1].chr))
				ref_len += H[nH-1].e;
			H[nH].chr = copy_string(nam);
			H[nH].b = i;
			H[nH].e = j;
			// assumes last event per chrom. is a het. interval
			het_len += (j - i);
			++nH;
		}
		fclose(fp);
	}
	ref_len += H[nH-1].e;

	// populations must be disjoint
	saw[0] = saw[1] = saw[2] = saw[3] = 0;
	for (i = 9; i < argc; ++i) {
		if (sscanf(argv[i], "%d:%d:%s", &j, &k, nam) != 3)
			fatalf("not like 13:2:fred : %s", argv[i]);
		if (k < 0 || k > 3)
			fatalf("not population 0, 1, 2 or 3: %s", argv[i]);
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
		} else {	// in a source population
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
	pop3 = saw[3];

	// start the output file of text
	for (k = 1; k <= 3; ++k) {
		if (k == 3 && !pop3)
			break;
		fprintf(out, "source population %d (state %d):", k, k);
		for (i = 0; i < nI; ++i)
			if (C[i].pop == k)
				fprintf(out, " %s", C[i].name);
		fprintf(out, "\n\n");
	}
	if (pop3) {
		fprintf(out, "state 4 is heterozygous for populations 1 and 2\n");
		fprintf(out,
		   "state 5 is heterozygous for populations 1 and 3\n");
		fprintf(out,
		   "state 6 is heterozygous for populations 2 and 3\n");
	} else
		fprintf(out, "state 3 is heterozygous for populations 1 and 2\n");
	fprintf(out, "\nswitch penalty = %2.2f\n", switch_penalty);
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

	if (ref_len)
		fprintf(out,
		  "%lld of %lld reference bp (%1.1f%%) are heterochromatin\n\n",
		  het_len, ref_len, 100.0*(float)het_len/(float)ref_len);

	// write fractions in each state to the output text file
	for (j = 0; j < 7; ++j)
		tot[j] = 0.0;
	fprintf(out, "individual:");
	fprintf(out, "\tstate 1\tstate 2\tstate 3");
	if (pop3)
		fprintf(out, "\tstate 4\tstate 5\tstate 6");
	fprintf(out, "\t  pop 1\t  pop 2");
	if (pop3)
		fprintf(out, "\t  pop 3");
	putc('\n', out);
	for (i = 0; i < nG; ++i) {
		N = A[i].x[1] + A[i].x[2] + A[i].x[4];
		if (pop3)
			N += A[i].x[3] + A[i].x[5] + A[i].x[6];
		N /= 100.0;
		fprintf(out, "%s:", A[i].name);
		if (strlen(A[i].name) < 7)
			putc('\t', out);
		for (j = 1; j < 7; ++j)
			if (pop3 || j == 1 || j == 2 || j == 4) {
				tot[j] += (keep[j] = A[i].x[j]);
				fprintf(out, "\t %5.1f%%", keep[j]/N);
			}
		keep[1] += 0.5*keep[4];
		keep[2] += 0.5*keep[4];
		if (pop3) {
			keep[1] += 0.5*keep[5];
			keep[2] += 0.5*keep[6];
			keep[3] += 0.5*(keep[5]+keep[6]);
		}

		fprintf(out, "\t %5.1f%%\t %5.1f%%", keep[1]/N, keep[2]/N);
		if (pop3)
			fprintf(out, "\t %5.1f%%", keep[3]/N);
			
		putc('\n', out);
	}
	if (nG > 1) {
		fprintf(out, "\naverage: ");
		N = tot[1] + tot[2] + tot[4];
		if (pop3)
			N += (tot[3] + tot[5] + tot[6]);
		N /= 100.0;
		for (j = 1; j < 7; ++j) {
			if (pop3 || j == 1 || j == 2 || j == 4)
				fprintf(out, "\t %5.1f%%", tot[j]/N);
		}
		xx = tot[1] + 0.5*tot[4];
		yy = tot[2] + 0.5*tot[4];
		if (pop3) {
			xx += 0.5*tot[5];
			yy += 0.5*tot[6];
		}
		fprintf(out, "\t %5.1f%%\t %5.1f%%", xx/N, yy/N);
		if (pop3)
			fprintf(out, "\t %5.1f%%",
			  (tot[3] + 0.5*tot[5] + 0.5*tot[6])/N);
		putc('\n', out);
	}

	return 0;
}
