#include "lib.h"

#define MAX_EVAL 1000

float E[MAX_EVAL];
int nE;

int main (int argc, char **argv) {
	FILE *fp;
	char buf[500];
	int i;
	float tot;

	fp = (argc== 1 ? stdin : ckopen(argv[1], "r"));
	while (fgets(buf, 500, fp)) {
		if (nE >= MAX_EVAL)
			fatal("Too many eigenvalues");
		E[nE++] = atof(buf);
	}
	for (tot = 0.0, i = 0; i < nE; ++i)
		tot += E[i];
	printf("Percentage explained by eigenvectors:\n");
	for (i = 0 ; i < nE && E[i] > 0.0; ++i)
		printf("%d: %1.1f%%\n", i+1, 100.0*(float)E[i]/tot);
	return 0;
}
