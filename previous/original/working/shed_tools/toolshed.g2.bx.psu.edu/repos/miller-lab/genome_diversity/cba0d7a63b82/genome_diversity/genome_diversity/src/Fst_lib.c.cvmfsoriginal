// This file contains three procedures for computing different variants of Fst.

#include "Fst_lib.h"

/* For each of wright, weir and reich, args 1 and 2 are counts for one
*  allele in the two populations; args 3 and 4 are counts for the other allele.
*  The numerator and denominator of the computed Fst are returned through
*  args 5 and 6.
*/

void wright(int A1, int A2, int B1, int B2, double *N, double *D) {
	double a1 = A1, a2 = A2, b1 = B1, b2 = B2, n1, n2, p1, p2;

	double
	  p, // frequency in the pooled population
	  H_ave, // average of HWE heterogosity in the two populations
	  H_all; // HWE heterozygosity in the pooled popuations

	n1 = a1+b1;
	n2 = a2+b2;
	if (n1 == 0.0 || n2 == 0.0) {
		// let the calling program handle it
		*N = *D = 0.0;
		return;
	}
	p1 = a1/n1;
	p2 = a2/n2;
	H_ave = p1*(1.0 - p1) + p2*(1.0 - p2);
	p = (p1 + p2)/2.0;
	H_all = 2.0*p*(1.0 - p);
	*N = H_all - H_ave;
	*D = H_all;
}

void weir(int A1, int A2, int B1, int B2, double *N, double *D) {
	double a1 = A1, a2 = A2, b1 = B1, b2 = B2, n1, n2, p1, p2,
	  n_tot, p_bar, nc, MSP, MSG;

	n1 = a1+b1;
	n2 = a2+b2;
	if (n1 == 0.0 || n2 == 0.0) {
		// let the calling program handle it
		*N = *D = 0.0;
		return;
	}
	n_tot = n1 + n2; 
	p1 = a1/n1;
	p2 = a2/n2;

	MSG = (n1*p1*(1.0-p1) + n2*p2*(1.0-p2))/(n_tot-2.0);
	p_bar = (n1*p1 + n2*p2)/n_tot;
	MSP = n1*(p1-p_bar)*(p1-p_bar) + n2*(p2-p_bar)*(p2-p_bar);
        nc = n_tot - (n1*n1 + n2*n2)/n_tot;
	*N = MSP - MSG;
	*D = MSP + (nc-1)*MSG;
}

void reich(int A1, int A2, int B1, int B2, double *N, double *D) {
	double a1 = A1, a2 = A2, b1 = B1, b2 = B2, n1, n2,  h1, h2, x;

	n1 = a1+b1;
	n2 = a2+b2;
	if (n1<=1 || n2<=1) {
		// let the calling program handle it
		*N = *D = 0.0;
		return;
	}
	h1 = (a1*(n1-a1)) / (n1*(n1-1));
	h2 = (a2*(n2-a2)) / (n2*(n2-1));
	x = a1/n1 - a2/n2;
	*N = x*x - h1/n1 - h2/n2;
	*D = *N + h1 + h2;
}
