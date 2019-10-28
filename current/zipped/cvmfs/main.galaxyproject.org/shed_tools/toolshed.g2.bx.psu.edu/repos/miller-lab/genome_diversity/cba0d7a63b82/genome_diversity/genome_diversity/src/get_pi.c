/* get_pi -- compute piNon, piSyn, thetaNon and thetaSyn
*
*  argv[1] -- SAPs file
*  argv[2] -- SNPs file
*  argv[3] -- covered intervals file
*  argv[4], argv[5], ... -- starting columns in the SNP file for the chosen
*    individuals
*
*  We assume that lines of argv[1], argv[2] and argv[3] start with a scaffold
*  name and a scaffold position, and that they are sorted on those two fields.
*  The 4th entry in an interval line gives the reference chromosome. We ignore
*  unnumbered chromosome, e.g., "chrX".
*
*  Output:
*    the number of nonsyn SNPs, total number of nonsynon sites,  piNon,
*    the number of synon SNPs, total number of synon sites, piSyn, plus
*    total length of covered intervals, thetaNon, thetaSyn
*
* What it does on Galaxy
The tool computes values that estimate some basic parameters
*/  

#include "lib.h"

// END_FILE should be larger than any real scaffold number
#define END_FILE 1000000000 // scaffold "number" signifying end-of-file
#define BUF_SIZE 50000	// maximum length of a SNP-file row

int col[10000], nC; // columns containing the chosen genotypes

// scaffold numbers and positions of the current SAP, SNP, and interval
int nbr_SAP, nbr_SNP, nbr_interv, pos_SAP, pos_SNP, beg, end, columns, debug;

// current SNP row, the variant amino acids of the SAP, interval's reference chr
char snp[BUF_SIZE], A[100], B[100], chr[100];

// number of nonsynon and snynon sites in the current interval
float all_non, all_syn;

// return the number of chromosome pairs that differ at a SNP
int diff_pair() {
	int i, j, X[1000];
	char *p, *z = "\t\n";

	// set X[i] = atoi(i-th word of SNP-file row), base 1
	for (i = 1, p = strtok(snp, z); p != NULL;
	  ++i, p = strtok(NULL, z))
		X[i] = atoi(p);
	// add genotypes to count the reference allele
	for (j = i = 0; i < nC; ++i)
		j += X[col[i]];
	// of the 2*nC chromosomes, j have the reference nucleotide
	if (debug)
		printf("get_pi: j = %d, return %d\n", j, j*(2*nC-j));
	return j*(2*nC-j);
}

// translate e.g. "scaffold234" to the integer 234
int name2nbr(char *s) {
	if (same_string(s, "chrX"))
		return 1000;
	if (same_string(s, "chrY"))
		return 1001;
	while (isalpha(*s))
		++s;
	return atoi(s);
}

// does one scaffold-position pair strictly precede another
int before(int nbra, int posa, int nbrb, int posb) {
	return (nbra < nbrb || (nbra == nbrb && posa < posb));
}

// get a SAP; set A and B; set nbr_SAP = END_FILE for end-of-file
void get_SAP(FILE *fp) {
	char buf[500], scaf_name[100];
	int old_nbr = nbr_SAP, old_pos = pos_SAP;

	if (nbr_SAP >= END_FILE)
		return;
	if (!fgets(buf, 500, fp)) {
		nbr_SAP = END_FILE;
		return;
	}
	while (buf[0] == '#')
		if (!fgets(buf, 500, fp)) {
			nbr_SAP = END_FILE;
			return;
		}
	if (columns == 8) {
		if (sscanf(buf, "%s %d %*s %*s %*s %*s %s %*s %s",
		    scaf_name, &pos_SAP, A, B) != 4)
			fatalf("bad SAP : %s", buf);
	} else if (columns == 5) {
		if (sscanf(buf, "%s %d %*s %*s %s %*s %s",
		    scaf_name, &pos_SAP, A, B) != 4)
			fatalf("bad SAP : %s", buf);
	} else
		fatalf("get_SAP: columns = %d", columns);
	nbr_SAP = name2nbr(scaf_name);
	if (before(nbr_SAP, pos_SAP, old_nbr, old_pos))
		fatalf("SAP at scaf%d %d is out of order", nbr_SAP, pos_SAP);
	if (debug)
		printf("SAP: scaf%d %d\n", nbr_SAP, pos_SAP);
}

// get a SNP
void get_SNP(FILE *fp) {
	char scaf_name[100];
	int old_nbr = nbr_SNP, old_pos = pos_SNP;

	if (nbr_SNP >= END_FILE)
		return;
	if (!fgets(snp, BUF_SIZE, fp)) {
		nbr_SNP = END_FILE+1;
		return;
	}
	while (snp[0] == '#')
		if (!fgets(snp, 500, fp)) {
			nbr_SNP = END_FILE+1;
			return;
		}
	if (sscanf(snp, "%s %d", scaf_name, &pos_SNP) != 2)
		fatalf("bad SNP : %s", snp);
	nbr_SNP = name2nbr(scaf_name);
	if (before(nbr_SNP, pos_SNP, old_nbr, old_pos)) {
		fprintf(stderr, "seq%d %d before seq%d %d\n",
		  nbr_SNP, pos_SNP, old_nbr, old_pos);
		fatalf("SNP at sequence %d %d is out of order", nbr_SNP, pos_SNP);
	}
	if (debug)
		printf("SNP: scaf%d %d\n", nbr_SNP, pos_SNP);
}

// expand fractions .333 and .666 to full double-precision accuracy
double grow(float x) {
	int chop = x;
	float remain;
	double d, third = (double)1/(double)3;

	d = (double)chop;
	remain = x - (float)chop;
	if (0.1 < remain)
		d += third;
	if (0.5 < remain)
		d += third;
	return d;
}

// read an interval; update tot_non and tot_syn
int get_interval(FILE *fp) {
	char buf[500], scaf_name[100], tmp[500], *t, *z = " \t\n";
	int old_nbr = nbr_interv, old_end = end;

	if (!fgets(buf, 500, fp))
		return 0;
	while (buf[0] == '#')
		if (!fgets(buf, 500, fp))
			return 0;
	if (columns == 0) {
		strcpy(tmp, buf);
		for (columns = 0, t = strtok(tmp, z); t != NULL;
		     ++columns, t = strtok(NULL, z))
			;
	} 
	if (columns != 5 && columns != 8)
		fatalf("interval table has %d columns", columns);
	if (columns == 8 && sscanf(buf, "%s %d %d %s %*s %*s %f %f",
	    scaf_name, &beg, &end, chr, &all_non, &all_syn) != 6)
		fatalf("bad interval : %s", buf);
	if (columns == 5) {
		if (sscanf(buf, "%s %d %d %f %f",
		    scaf_name, &beg, &end, &all_non, &all_syn) != 5)
			fatalf("bad interval : %s", buf);
		strcpy(chr, scaf_name);
	}
	nbr_interv = name2nbr(scaf_name);
	if (before(nbr_interv, beg, old_nbr, old_end))
		fatalf("interval at %s %d is out of order", scaf_name, beg);
	if (debug)
		printf("int: scaf%d %d-%d\n", nbr_interv, beg, end);
		
	return 1;
}

int main(int argc, char **argv) {
	FILE *fp1, *fp2, *fp3;
	int i, nint, nsap, no_sap, no_snp, no_chr, nsyn, nnon, d, tot_len;
	float non, syn, x;
	double tot_non = 0.0, tot_syn = 0.0,	// total sites in the intervals
	  factor;

	// process command-line arguments
	if (same_string(argv[argc-1], "debug")) {
		debug = 1;
		--argc;
	}
	if (argc < 5)
		fatal("args: SAPs SNPs intervals individual1 ... [debug]");
	fp1 = ckopen(argv[1], "r");
	fp2 = ckopen(argv[2], "r");
	fp3 = ckopen(argv[3], "r");
	for (i = 4; i < argc; ++i)
		col[i-4] = atoi(argv[i]) + 2;
	nC = argc - 4;

	// loop over the intervals
	tot_len = no_sap = nsap = no_snp = no_chr = nsyn = nnon = 0;
	non = syn = 0.0;
	for (nint = 0; get_interval(fp3); ++nint) {
		if (strncmp(chr, "chr", 3) == 0 && !isdigit(chr[3])) {
			++no_chr;
			continue;
		}
		tot_len += (end - beg);
		// expand e.g. .333 to .3333333..
		tot_non += grow(all_non);
		tot_syn += grow(all_syn);

		// skip SAPs coming before this interval
		while (before(nbr_SAP, pos_SAP, nbr_interv, beg))
			get_SAP(fp1);
		// loop over SAPs in this inteval
		while (before(nbr_SAP, pos_SAP, nbr_interv, end)) {
			++nsap;

			// look for this SAP in the SNP file
			while (before(nbr_SNP, pos_SNP, nbr_SAP, pos_SAP)) {
				if (nbr_SNP == nbr_interv && pos_SNP >= beg)
					++no_sap;
				get_SNP(fp2);
			}

			// is this the SAP we looked for?
			if (nbr_SAP == nbr_SNP && pos_SAP == pos_SNP) {
				d = diff_pair();
				if (A[0] == B[0]) {
					++nsyn;
					syn += (float)d;
				} else {
					++nnon;
					non += (float)d;
				}
				get_SNP(fp2);
			} else
				++no_snp;
			get_SAP(fp1);
		}
		// process SNPs in the interval but not in the SAP file
		while (before(nbr_SNP, pos_SNP, nbr_interv, end)) {
			if (nbr_SNP == nbr_interv && pos_SNP >= beg)
				++no_sap;
			get_SNP(fp2);
		}
	}

	// there are x = (2*nC choose 2) pairs of chromosomes
	x = (float)(nC*(2*nC-1));
	non /= x;
	syn /= x;
	printf("%d intervals\n", nint);
	if (no_chr)
		printf("ignored %d interval%s on unnumbered chromosomes, like chrX\n",
		  no_chr, no_chr == 1 ? "" : "s");
	printf("%d SNPs, %d nonsyn, %d synon\n", nsap, nnon, nsyn);
	if (no_sap)
		printf("%d SNPs in an interval are not in the SAP table\n",
		  no_sap);
	if (no_snp)
		printf("%d SNPs in an interval are not in the SNP table\n",
		  no_snp);
	printf("nonsynon: %4.3f/%4.3f = %6.5f%%\n",
	  non, tot_non, 100.0*non/tot_non);
	printf("synon: %4.3f/%4.3f = %6.5f%%\n",
	  syn, tot_syn, 100.0*syn/tot_syn);
	for (factor = 0.0, i = 1; i < 2*nC; ++i)
		factor += (1.0/(double)i);
	factor *= (double)tot_len/100.0;
	printf("%d covered bp, thetaNon = %6.5f%%, thetaSyn = %6.5f%%\n",
	 tot_len, (double)nnon/factor, (double)nsyn/factor);
	return 0;
}
