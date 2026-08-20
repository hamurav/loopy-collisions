CC     ?= cc
CFLAGS ?= -O3 -pthread
# On x86 add -march=native for ~10-15%.  Do NOT use -fopenmp: this is pthreads.

all: loopy_collisions comp_check gen_graphs mg_collisions gen_multi

loopy_collisions: loopy_collisions.c ; $(CC) $(CFLAGS) -o $@ $<
comp_check:       comp_check.c       ; $(CC) $(CFLAGS) -o $@ $<
gen_graphs:       gen_graphs.c       ; $(CC) -O3 -o $@ $<
mg_collisions:    mg_collisions.c    ; $(CC) $(CFLAGS) -o $@ $<
gen_multi:        gen_multi.c        ; $(CC) -O3 -o $@ $<

# quick self-test: reproduces the published n <= 9 answers in a few seconds
check: all
	./gen_graphs 8 > /tmp/c8.g6 && ./loopy_collisions -q < /tmp/c8.g6 | grep -E 'SUMMARY|VERDICT'
	./gen_graphs 9 > /tmp/c9.g6 && ./loopy_collisions -q < /tmp/c9.g6 | grep -E 'SUMMARY|VERDICT'
	python3 mn_example.py
	./gen_graphs 5 | ./gen_multi 3 | ./mg_collisions -q | grep -E 'SUMMARY|VERDICT'
	./gen_graphs 6 | python3 min_edges.py --maxm 9

clean: ; rm -f loopy_collisions comp_check gen_graphs mg_collisions gen_multi
.PHONY: all check clean
