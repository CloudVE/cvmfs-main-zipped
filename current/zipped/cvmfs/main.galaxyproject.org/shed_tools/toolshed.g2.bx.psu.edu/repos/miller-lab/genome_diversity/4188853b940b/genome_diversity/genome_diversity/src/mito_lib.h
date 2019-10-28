// mito_data.h -- header file for shared procedures to read SNPs and intervals
// for mitogenomes

#define MAX_SNP 20000
#define MAX_SAMPLE 100
struct snp {
	int pos;
	char *g; // genotypes - one character per specified mitogenome
} S[MAX_SNP], I[MAX_SNP];

// intervals associated with each specified mitogenome
struct interv {
	int b, e;
	struct interv *next;
};
int nM, min_cov, debug;

// mitogenomes
struct mito {
	char *name;
	int col;	// first column in the SNP table
	struct interv *intervals, *last;
} M[MAX_SAMPLE];

// get the adequately covered intervals for each specified individual;
// merge adjacent intervals
void get_intervals(char *filename);

// get the SNPs; for each SNP set the array of (first characters from)
// genotypes of the specified samples (individuals)
int get_variants(char *filename, struct snp *S, int refcol);
