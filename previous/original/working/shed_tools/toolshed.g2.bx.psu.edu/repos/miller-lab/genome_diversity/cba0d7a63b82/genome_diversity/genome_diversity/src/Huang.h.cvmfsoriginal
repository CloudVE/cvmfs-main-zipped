/* Find intervals of highest total score, i.e., such that adding postions to
*  either end will decrease the total. We use the method of Fig. 6 of the paper:
*  Xiaoqiu Huang, Pavel Pevzner, Webb Miller (1994) Parametric recomputing in
*  alignment graphs. Combinatorial Pattern Matching (Springer Lecture Notes in
*  Computer Science, 807), 87-101.
*
*  The input scores are in x[0], x[1], ..., x[n-1], but the output regions
*  are in R[1], R[2], ..., R[top]. R[i].Score is the total score of the i-th
*  (in order of position) positive-scoring interval of x, which consists of of
*  x[R[i].Lpos + 1] to x[R[i].Rpos].
*/
#define MAX_R 5000000

struct region {	// a consecutive (relative to the reference) run of SNPs
	double Lscore, Rscore, Score;
	int Lpos, Rpos, Lower;
} R[MAX_R];
int top;

void Huang(double *x, int n);
