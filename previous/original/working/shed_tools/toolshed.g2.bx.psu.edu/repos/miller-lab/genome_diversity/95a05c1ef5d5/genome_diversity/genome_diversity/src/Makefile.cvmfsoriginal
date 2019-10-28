CC = gcc
COPT = -O2
CWARN = -W -Wall
CFLAGS = $(COPT) $(CWARN)
INSTALL_DIR = ../bin

TARGETS = admix_prep aggregate coords2admix coverage dist_mat dpmix \
          eval2pct filter_snps Fst_ave Fst_column get_pi sweep

all: $(TARGETS)

install: $(TARGETS)
	if [ ! -d "$(INSTALL_DIR)" ]; then mkdir -p "$(INSTALL_DIR)"; fi
	cp $(TARGETS) $(INSTALL_DIR)

admix_prep: admix_prep.c lib.c
	$(CC) $(CFLAGS) $^ -o $@

aggregate: aggregate.c lib.c
	$(CC) $(CFLAGS) $^ -o $@

coords2admix: coords2admix.c lib.c
	$(CC) $(CFLAGS) $^ -o $@

coverage: coverage.c lib.c
	$(CC) $(CFLAGS) $^ -o $@

dist_mat: dist_mat.c lib.c
	$(CC) $(CFLAGS) $^ -o $@

dpmix: dpmix.c lib.c
	$(CC) $(CFLAGS) $^ -o $@

eval2pct: eval2pct.c lib.c
	$(CC) $(CFLAGS) $^ -o $@

filter_snps: filter_snps.c lib.c
	$(CC) $(CFLAGS) $^ -o $@

Fst_ave: Fst_ave.c Fst_lib.c lib.c
	$(CC) $(CFLAGS) $^ -o $@

Fst_column: Fst_column.c Fst_lib.c lib.c
	$(CC) $(CFLAGS) $^ -o $@

get_pi: get_pi.c lib.c
	$(CC) $(CFLAGS) $^ -o $@

sweep: sweep.c lib.c Huang.c
	$(CC) $(CFLAGS) $^ -o $@

.PHONY: clean

clean:
	rm -f $(TARGETS)
