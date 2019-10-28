// Find highest scoring intervals, as discussed in Huang.h.

#include "lib.h"
#include "Huang.h"

void Huang(double x[], int n) {
	double Score, oldScore;
	int v, L, i;

	top = 0;	// don't use location 0, so as to follow Fig. 6
	for (Score = 0.0, v = 0; v < n; ++v) {
		oldScore = Score;
		Score += x[v];
		if (x[v] < 0)
			continue;
		if (top > 0 && R[top].Rpos == v-1) {
			// add edge to top subpath
			R[top].Rpos = v;
			R[top].Rscore = Score;
		} else {
			// create a one-edge subpath
			++top;
			if (top >= MAX_R)
				fatal("In Haung(), top is too big");
			R[top].Lpos = v-1;
			R[top].Lscore = oldScore;
			R[top].Rpos = v;
			R[top].Rscore = Score;
			R[top].Lower = top-1;
			while ((L = R[top].Lower) > 0 &&
			  R[L].Lscore > R[top].Lscore)
				R[top].Lower = R[L].Lower;
		}
		// merge subpaths
		while (top > 1 && (L = R[top].Lower) > 0 &&
		    R[L].Rscore <= R[top].Rscore) {
			R[L].Rpos = R[top].Rpos;
			R[L].Rscore = R[top].Rscore;
			top = L;
		}
	}
	for (i = 1; i <= top; ++i)
		R[i].Score = R[i].Rscore - R[i].Lscore;
}
