CC = gcc
COPT = -O2
CWARN = -W -Wall
CFLAGS = $(COPT) $(CWARN)
INSTALL_DIR = ../bin

TARGETS = admix_prep coords2admix coverage dist_mat dpmix eval2pct \
          Fst_ave Fst_column pop sweep

all: $(TARGETS)

install: $(TARGETS)
	if [ ! -d "$(INSTALL_DIR)" ]; then mkdir -p "$(INSTALL_DIR)"; fi
	cp $(TARGETS) $(INSTALL_DIR)

admix_prep: admix_prep.c lib.c
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

Fst_ave: Fst_ave.c Fst_lib.c lib.c
	$(CC) $(CFLAGS) $^ -o $@

Fst_column: Fst_column.c Fst_lib.c lib.c
	$(CC) $(CFLAGS) $^ -o $@

pop: pop.c lib.c
	$(CC) $(CFLAGS) $^ -o $@

sweep: sweep.c lib.c Huang.c
	$(CC) $(CFLAGS) $^ -o $@

.PHONY: clean

clean:
	rm -f $(TARGETS)
