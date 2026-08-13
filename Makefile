CC     ?= cc
CFLAGS ?= -O3 -pthread
# On x86 add -march=native for ~10-15%.  Do NOT use -fopenmp: this is pthreads.

all: loopy_collisions comp_check gen_graphs

loopy_collisions: loopy_collisions.c ; $(CC) $(CFLAGS) -o $@ $<
comp_check:       comp_check.c       ; $(CC) $(CFLAGS) -o $@ $<
gen_graphs:       gen_graphs.c       ; $(CC) -O3 -o $@ $<

# quick self-test: reproduces the published n <= 9 answers in a few seconds
check: all
	./gen_graphs 8 > /tmp/c8.g6 && ./loopy_collisions -q < /tmp/c8.g6 | grep -E 'SUMMARY|VERDICT'
	./gen_graphs 9 > /tmp/c9.g6 && ./loopy_collisions -q < /tmp/c9.g6 | grep -E 'SUMMARY|VERDICT'

clean: ; rm -f loopy_collisions comp_check gen_graphs
.PHONY: all check clean
