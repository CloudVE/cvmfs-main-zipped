// coords2admix -- add projections onto chords to information about
// coordinates in PCA plots

#include "lib.h"

#define MAX_POP 1000
struct pop {
	char *name;
	float x, y;
} P[MAX_POP];
int nP;

int main(int argc, char **argv) {
	FILE *fp;
	char buf[500], x[100], y[100], z[100], cur_pop[100];
	int ncur, i, j, k;
	float eig1, eig2, tot_x = 0.0, tot_y = 0.0, x1, y1, x2, y2, a, b, c, d;

	if (argc == 1)
		fp = stdin;
	else if (argc == 2)
		fp = ckopen(argv[1], "r");
	else
		fatal("optional arg: smartpca coordinates");

	if (!fgets(buf, 500, fp))
		fatal("empty set of coordinates");
	if (sscanf(buf, "%s %s %s", x, y, z) != 3 ||
	    !same_string(x, "#eigvals:"))
		fatalf("cannot find eigenvalues: %s", buf);
	printf("%s", buf);
	eig1 = atof(y);
	eig2 = atof(z);
	//printf("eig1 = %f, eig2 = %f\n", eig1, eig2);
	
	strcpy(cur_pop, "");
	ncur = 0;
	while (fgets(buf, 500, fp)) {
		if (sscanf(buf, "%*s %s %s %s", x, y, z) != 3)
			fatalf("gag: %s", buf);
		printf("%s", buf);
		if (!same_string(cur_pop, z)) {
			if (ncur > 0) {
				P[nP].name = copy_string(cur_pop);
				P[nP].x = tot_x/ncur;
				P[nP].y = tot_y/ncur;
				++nP;
			}
			ncur = 1;
			strcpy(cur_pop, z);
			tot_x = atof(x);
			tot_y = atof(y);
		} else {
			++ncur;
			tot_x += atof(x);
			tot_y += atof(y);
		}
	}
	P[nP].name = copy_string(cur_pop);
	P[nP].x = tot_x/ncur;
	P[nP].y = tot_y/ncur;
	++nP;

/*
for (i = 0; i < nP; ++i)
printf("%s %f %f\n", P[i].name, P[i].x, P[i].y);
*/

	// loop over pairs of populations
	for (i = 0; i < nP; ++i) {
		x1 = eig1*P[i].x;
		y1 = eig2*P[i].y;
		for (j = i+1; j < nP; ++j) {
			printf("\nprojection along chord %s -> %s\n",
			  P[i].name, P[j].name);
			x2 = eig1*P[j].x;
			y2 = eig2*P[j].y;
			c = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
			for (k = 0; k < nP; ++k)
				if (k != i && k != j) {
					a = eig1*P[k].x;
					b = eig2*P[k].y;
					d = (x2-x1)*(a-x1) + (y2-y1)*(b-y1);
					printf("  %s: %f\n", P[k].name, d/c);
				}
		}
	}

	return 0;
}

