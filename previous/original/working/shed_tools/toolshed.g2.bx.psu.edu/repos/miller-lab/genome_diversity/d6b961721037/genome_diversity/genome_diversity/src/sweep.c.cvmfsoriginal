/* sweep -- find regions of the genome with high scores (e.g., Fst scores).
*
*  argv[1] -- file containing a Galaxy table
*  argv[2] -- column number (base-1) for the chromosome name
*  argv[3] -- column number for the (base-0) chromosomal position
*  argv[4] -- column number for a score for the position
*  argv[5] -- a percentage, such as "95", or a raw score, such as "=0.9".
*  argv[6] -- the number of randomizations (shuffles) of the scores
*  argv[7] -- [optional] if present and non-zero, report SNPs
*
*  The program first determines a threshold such that the stated percentage
*  of the scores are below that threshold (or uses the provided number if
*  argv[5] starts with "=").  The program subtracts the threshold
*  from each score, then looks for maximal-scoring runs of SNPs, i.e., where
*  adding or subtracting SNPs from an end of then run always decreases the
*  total score. These regions are printed in order of descreasing total score.
*  To determine a cutoff for the printed regions, the programs takes the maximum
*  score over all regions observed in a specified number of shuffles of the
*  list of scores. If argv[6] = 0, then all maximal-scoring runs of at least
*  4 table entries are printed.

What it does on Galaxy
The user selects a SNP table and specifies the columns containing (1) chromosome, (2) position, (3) scores (such as an Fst-value for the SNP), (4) a percentage or raw score for the "cutoff" and (5) the number of times the data should be radomized (only intervals with score exceeding the maximum for the randomized data are reported). If a percentage (e.g. 95%) is specified for #3, then that percentile of the scores is used as the cutoff; this may not work well if many SNPs have the same score. The program subtracts the cutoff from every score, then finds genomic intervals (i.e., consecutive runs of SNPs) whose total score cannot be increased by adding or subtracting one or more SNPs at the ends of the interval.
*/

#include "lib.h"
#include "Huang.h"

// maximum number of rows in any processed table
#define MANY 20000000
#define BUF_SIZE 5000
#define MAX_WINDOW 1000000

double X[MANY];	// holds all scores
int nX;

// position-score pairs for a single chromosome
struct score {
	int pos;
	double x; // original score, then shifted score
} S[MANY];
int nS;

struct snp {
	int pos;
	double x;
	struct snp *next;
};

// structure to hold the maximum-scoring chromosomal intervals
struct sweep {
	float score;
	char *chr;
	int b, e;
	struct snp *snps;
} W[MAX_WINDOW];
int nW;

// return the linked list of SNPs in positions b to e
struct snp *add_snps(int b, int e) {
	struct snp *first = NULL, *last = NULL, *new;
	int i;
	for (i = b; i <= e; ++i)
		if (S[i].pos >= 0) {
			new = ckalloc(sizeof(*new));
			new->pos = S[i].pos;
			new->x = S[i].x;
			new->next = NULL;
			if (first == NULL)
				first = new;
			else
				last->next = new;
			last = new;
		}
	return first;
}

// given a table row, return a pointer to the item in a particular column
char *get_col(char *buf, int col) {
	static char temp[BUF_SIZE], *p;
	int i;
	char *z = " \t\n";

	strcpy(temp, buf);
	for (p = strtok(temp, z), i = 1; *p && i < col;
	     p = strtok(NULL, z), ++i)
		;
	if (p == NULL)
		fatalf("no column %d in %s", col, buf);
	return p;
}

// fill S[] with position-score pairs for the next chromosome
// return 0 for EOF
int get_chr(FILE *fp, int chr_col, int pos_col, int score_col, char *chr) {
	static char buf[BUF_SIZE];
	static int init = 1;
	char *status;

	if (init) {
		while ((status = fgets(buf, BUF_SIZE, fp)) != NULL &&
		  buf[0] == '#')
			;
		if (status == NULL)
			fatal("empty table");
		init = 0;
	}
	if (buf[0] == '\0')
		return 0;
	
	if (buf[0] == '#')
		fatal("cannot happen");
	strcpy(chr, get_col(buf, chr_col));
	S[0].pos = atoi(get_col(buf, pos_col));
	S[0].x = atof(get_col(buf, score_col));
	for (nS = 1; ; ++nS) {
		if (!fgets(buf, BUF_SIZE, fp)) {
			buf[0] = '\0';
			return 1;
		}
		if (!same_string(chr, get_col(buf, chr_col)))
			break;
		S[nS].pos = atoi(get_col(buf, pos_col));
		S[nS].x = atof(get_col(buf, score_col));
	}
	return 1;
}

// for sorting genomic intervals by *decreasing* score
int Wcompar(struct sweep *a, struct sweep *b) {
	float y = a->score, z = b->score;

	if (y > z)
		return -1;
	if (y < z)
		return 1;
	return 0;
}

// for sorting an array of scores into increasing order
int fcompar(double *a, double *b) {
	if (*a < *b)
		return -1;
	if (*a > *b)
		return 1;
	return 0;
}

/* shuffle the values S[0], S[1], ... , S[nscores-1];
*  Uses Algorithm P in page 125 of "The Art of Computer Programming (Vol II)
*  Seminumerical Programming", by Donald Knuth, Addison-Wesley, 1971.
*/
void shuffle_scores() {
	int i, j;
	double temp;

	for (i = nX-1; i > 0; --i) {
		// swap what's in location i with location j, where 0 <= j <= i
		j = random() % (i+1);
		temp = X[i];
		X[i] = X[j];
		X[j] = temp;
	}
}

// return the best interval score (R[i] is the struct operated by Huang())
double best() {
	int i;
	double bestScore;

	Huang(X, nX);

	for (bestScore = 0.0, i = 1; i <= top; ++i) 
		bestScore = MAX(R[i].Score, bestScore);
	return bestScore;
}

int main(int argc, char **argv) {
	FILE *fp;
	char buf[BUF_SIZE], chr[100], *a;
	double shift = 0.0, cutoff;
	int i, b, e, chr_col, pos_col, score_col, nshuffle, snps = 0;
	struct snp *s;

	if (argc != 7 && argc != 8)
		fatal("args: table chr_col pos_col score_col threhold randomizations [SNPs]");

	// process command-line arguments
	chr_col = atoi(argv[2]);
	pos_col = atoi(argv[3]);
	score_col = atoi(argv[4]);
	a = argv[5];
	fp = ckopen(argv[1], "r");
	if (argc == 8)
		snps = atoi(argv[7]);
	if (isdigit(a[0])) {
		for (nX = 0; nX < MANY && fgets(buf, BUF_SIZE, fp); ) {
			if (buf[0] == '#') 
				continue;
			X[nX++] = atof(get_col(buf, score_col));
		}
		if (nX == MANY)
			fatal("Too many rows");
		qsort((void *)X, (size_t)nX, sizeof(double),
		  (const void *)fcompar);
		shift = X[atoi(a)*nX/100];
		rewind(fp);
	} else if (a[0] == '=')
		shift = atof(a+1);

//fprintf(stderr, "shift = %4.3f\n", shift);
	nshuffle = atoi(argv[6]);
	if (nshuffle == 0)
		cutoff = 0;
	else {
		for (nX = 0; nX < MANY && fgets(buf, BUF_SIZE, fp); ) { 
			if (buf[0] == '#')
				continue;
			X[nX++] = atof(get_col(buf, score_col)) - shift;
		}
		if (nX == MANY)
			fatal("Too many rows");
		for (cutoff = 0.0, i = 0; i < nshuffle; ++i) {
			shuffle_scores();
			cutoff = MAX(cutoff, best());
		}
		rewind(fp);
	}
//fprintf(stderr, "cutoff = %4.3f\n", cutoff);

	// loop over chromosomes;
	// start by getting the chromosome's scores
	while (get_chr(fp, chr_col, pos_col, score_col, chr)) {
		// subtract shift from the scores
		for (i = 0; i < nS; ++i)
			X[i] = S[i].x - shift;

		// find the maximum=scoring regions
		Huang(X, nS);
	
		// save any regions with >= 4 points and score >= cutoff
		for (i = 0; i <= top; ++i) {
			if (nW >= MAX_WINDOW)
				fatalf("too many windows");

			// get indices of the first and last SNP in the interval
			b = R[i].Lpos + 1;
			e = R[i].Rpos;

			// remove unmapped SNP position from intervals' ends
			while (b < e && S[b].pos == -1)
				++b;
			while (e > b && S[e].pos == -1)
				--e;

			// record intervals
			if (e - b < 3 || R[i].Score < cutoff)
				continue;
			W[nW].score = R[i].Score;
			W[nW].chr = copy_string(chr);
			W[nW].b = S[b].pos;
			W[nW].e = S[e].pos+1;	// Ws are half-open
			if (snps)
				W[nW].snps = add_snps(b, e);
			++nW;
		}
	}

	// sort by decreasing score
	qsort((void *)W, (size_t)nW, sizeof(W[0]), (const void *)Wcompar);

	for (i = 0; i < nW; ++i) {
		printf("%s\t%d\t%d\t%4.4f\n", 
			W[i].chr, W[i].b, W[i].e, W[i].score);
		for (s = W[i].snps; s; s = s->next)
			printf(" %d %3.2f\n", s->pos, s->x);
	}
	return 0;
}
