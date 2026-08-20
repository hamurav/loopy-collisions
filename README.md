# loopy-collisions

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21929427.svg)](https://doi.org/10.5281/zenodo.21929427)

Code accompanying **A. Kirillov, G. Nenashev, B. Shapiro and A. Vaintrob,
*The Loopy Polynomial: from Tutte's Universal V-Function to Bizonotopal
Geometry*** (Sections 4.3 and 4.4).

Two searches, sharing one engine.

**Simple graphs.** Exhaustive search for pairs of non-isomorphic connected
simple graphs with equal **loopy polynomial** `L_G`, equal **refined loopy
polynomial** `Lhat_G`, or equal **Tutte symmetric function** `XB_G`, checking
that the three induce the same partition. They do, through **n = 11** vertices
(1 006 700 565 connected graphs).

**Loopless multigraphs.** The same search with multiplicities, over 48 229 871
multigraphs. Here the three do *not* agree: `L` and `Lhat` still induce the
same partition, but `XB` (equivalently the U-polynomial) is strictly coarser.
The smallest witness has five vertices and ten edges, and answers a question of
Merino and Noble; `python3 mn_example.py` verifies it from the definitions.

## Quick start

```sh
make            # builds all five programs
make check      # reproduces the published small answers, a few seconds
```

`make check` should print, for `n = 9`:

```
SUMMARY L n=9 graphs=152 classes=65 pairs=65 edges=[13,23]
VERDICT  L-classes == Lhat-classes : YES
VERDICT  Lhat-classes == XB-classes: YES
VERDICT  L-classes == XB-classes   : YES
```

and, for the multigraph half, the `n = 5`, `mu = 3` row of the table below
together with

```
VERDICT  Lhat-classes == XB-classes  : NO  <-- COUNTEREXAMPLE
```

On macOS use the default `cc` (clang). Do **not** pass `-fopenmp`; the code
uses pthreads. Thread count comes from `$NTHREADS`, else the online CPU count.

## Running a full search

Graph generation uses `geng` from [nauty](https://pallini.di.uniroma1.it/)
(McKay–Piperno). `gen_graphs.c` is a self-contained fallback for `n <= 9` that
self-validates against OEIS A000088/A001349.

Collisions occur only within a fixed edge count, so the work splits into
independent buckets:

```sh
for m in $(seq 10 55); do
  geng -q -c 11 ${m}:${m} | ./loopy_collisions -q > res11/m$m.txt
done
```

Large buckets can be split further with `--shard i/N` (sharding is on the
stage-0 key, so members of a class always land in the same shard). Note that
each shard still reads the whole bucket into memory — to reduce peak RSS,
split the `.g6` file itself instead.

Cost, measured: **1.5 min** for `n = 10`, **7.5 h** for `n = 11`, peak memory
about 5 GB (set by the largest bucket, `m = 27`, ~9.6e7 graphs).
`timing.sh RESDIR` recovers per-bucket wall times afterwards from file mtimes.

## Verifying a finished run

```sh
python3 verify_run.py res11 res10     # second argument optional
```

Checks: every bucket and shard produced output; classes recounted from the raw
`CLASS` lines and reconciled with the `SUMMARY` lines; the L / Lhat / XB
partitions identical class by class; the edge-count profile symmetric under
`m -> C(n,2) - m`; and every class of the `n-1` run reappearing at `m + (n-1)`
under `G -> K_1 \/ G`.

Two things that have caused confusion:

* `SUMMARY graphs=` counts **survivors of the sieve**, not graphs read.
* `SUMMARY pairs=` counts **unordered pairs**, so a class of size 4 contributes
  6, not 1.

## Producing the collision lists

A run leaves its classes scattered over one file per edge bucket (and per shard
of the sharded buckets), interleaved for all three invariants. To collect them
into one canonical file:

```sh
python3 make_lists.py res11 -o data/collisions_n11.txt
python3 make_lists.py res10 -o data/collisions_n10.txt
```

This deduplicates, checks that the L-, Lhat- and XB-classes agree (and exits
nonzero, loudly, if they do not), sorts by `(m, size, members)` and writes a
header with the class counts. Output is one class per line:

```
n  m  graph6 [graph6 ...]
```

with the class size given by the number of graph6 fields, so sizes 2, 3, 4 all
share one format. This is what `comp_check` and `lift.py` read.

## Other tools

* `exact_check.py res11/all.txt` — recomputes every reported class exactly over
  **Z**, removing any false positive from the modular fingerprints.
* `comp_check.c` — complements every class and checks it is still a class
  (`XB_Gbar` is an explicit invertible transform of `XB_G`).
* `lift.py` — applies `G -> K_1 \/ G`, which sends a class at `(n, m)` to one at
  `(n+1, m+n)` of the same size:
  `python3 lift.py < data/collisions_n10.txt` should be a subset of
  `data/collisions_n11.txt`.
* `gp.py` — slow, independent reference implementations straight from the
  defining expansions, used to validate the fast code on small graphs.
* `mn_example.py` — verifies the five-vertex multigraph pair of the paper's
  Merino--Noble theorem from the definitions, exactly over **Z**, sharing no
  code with anything else here. Prints the two displays quoted in the paper.
* `min_edges.py` — exhaustive search bounded by edge count rather than
  multiplicity; establishes that no counterexample has fewer than nine edges.
  See the multigraph section below.

## Results

| n | connected graphs | classes | sizes | edge range |
|---:|---:|---:|---|---:|
| ≤7 | 853 | 0 | — | — |
| 8 | 11 117 | 8 | 8×2 | 13–15 |
| 9 | 261 080 | 65 | 65×2 | 13–23 |
| 10 | 11 716 571 | 1 285 | 1281×2, 4×4 | 14–32 |
| 11 | 1 006 700 565 | 22 499 | 22491×2, 4×3, 4×4 | 13–42 |

In every case the L-, Lhat- and XB-classes coincide exactly. The complete
lists are in `data/`, one file per `n`, in the format described above.

## Multigraphs

`mg_collisions.c` runs the same search on loopless multigraphs. It reads the
text format of nauty's `multig`,

```
n m  v1 w1 mult1  v2 w2 mult2  ...
```

and computes the `L`, `Lhat`, `XB` and `U` fingerprints of every graph read;
`e_G(S)` simply counts edges with multiplicity, and nothing else in the
algorithm changes. There is no sieve — the cheap filters would each need
re-justifying for multigraphs, and at these vertex counts `3^n` is negligible.
Since `multig` emits pairwise non-isomorphic graphs, any class of size >= 2 is
a genuine collision and no isomorphism testing is needed.

```sh
geng -q -c 6 | multig -m4 -T | ./mg_collisions      # with nauty
./gen_graphs 6 | ./gen_multi 4 | ./mg_collisions    # without

sh run_multi.sh 6 4                                 # bucketed by edge count
```

`gen_multi.c` is a self-contained stand-in for `multig`: for each simple graph
it computes Aut(G) by brute force over the `n!` vertex permutations and keeps
one multiplicity vector per orbit. It is the bottleneck (`mu^|E|` assignments
per underlying graph), so prefer nauty when you have it.

`run_multi.sh n mu` splits the work by edge count with `mg_collisions -m M`,
one file per bucket, and reuses `all.mg` if it is already there. **Use it for
anything large**: the unbucketed `(6,4)` run needs about 3 GB resident and
thrashes, while the bucketed one peaks well under 1 GB. `check_partial.sh`
reports the must-hold verdicts on a run still in progress.

### Results

| n | mu | multigraphs | L-classes | Lhat-classes | XB-classes | U-classes |
|---:|---:|---:|---:|---:|---:|---:|
| 4 | 9 | 45 210 | 0 | 0 | 0 | 0 |
| 5 | 2 | 712 | 0 | 0 | 0 | 0 |
| 5 | 3 | 10 364 | 0 | 0 | 12 | 12 |
| 5 | 4 | 88 985 | 0 | 0 | 58 | 58 |
| 5 | 5 | 530 657 | 0 | 0 | 166 | 166 |
| 5 | 6 | 2 431 555 | 0 | 0 | 360 | 360 |
| 6 | 2 | 24 576 | 15 | 15 | 30 | 30 |
| 6 | 3 | 1 590 368 | 500 | 500 | 1 130 | 1 130 |
| 6 | 4 | 43 477 490 | 4 848 | 4 848 | 10 164 | 10 164 |
| 7 | 2 | 2 275 616 | 1 213 | 1 213 | 1 469 | 1 469 |

The rows nest — `(5,2) ... (5,5)` sit inside `(5,6)`, and `(6,2), (6,3)` inside
`(6,4)` — so the four maximal rows account for **48 229 871** multigraphs in
all. Throughout:

* `L`-classes and `Lhat`-classes are **identical**, on every graph tested;
* `XB`-classes and `U`-classes are identical (Sarmiento's equivalence);
* the `Lhat`-classes are **strictly finer** than the `XB`-classes. This is the
  counterexample: on loopless multigraphs `XB` collapses classes that the loopy
  polynomial keeps apart.

Every `L`-class found is a pair; two `XB`-classes at `(6,4)` have three
members. `data/multigraphs_n6_mu4.txt` holds all 4 848 `L`-classes and all
10 164 `XB`-classes of the largest run, produced by

```sh
python3 make_mg_lists.py res_6_4 -o data/multigraphs_n6_mu4.txt
```

which also re-derives the verdicts from the raw `CLASS` lines rather than
trusting the ones the C program printed, and exits nonzero if a must-hold one
fails.

**Validation.** With `mu = 1` the multigraph pipeline reproduces the published
simple-graph answers exactly (at `n = 8`: 8 classes, edges 13–15, all four
invariants agreeing) — worth rerunning after any change to `mg_collisions.c`,
since it exercises the same code paths on a known answer. The five-vertex
counterexample itself was confirmed exactly over **Z** by `mn_example.py` and,
independently, by `gp.py`.

### How few edges can a counterexample have?

The table above bounds the *multiplicity*. `min_edges.py` bounds the number of
**edges** instead, which is what settles minimality:

```sh
./gen_graphs 6 | python3 min_edges.py --maxm 9      # seconds
geng -q -c 8   | python3 min_edges.py --maxm 9      # with nauty
```

It reads connected simple graphs in graph6 as the possible *underlying* graphs,
enumerates every multiplicity vector with total at most `--maxm`, one per
`Aut(G)`-orbit, and computes `L` and `XB` exactly over **Z** from their
defining expansions — no fingerprints, no modular arithmetic, no code shared
with the C programs.

This answers the question completely for `m <= 9`. A connected graph with `m`
edges has at most `m+1` vertices, so a counterexample with at most eight edges
would live on at most nine vertices; and for `m = 9` the only case beyond
`n = 9` is a tree on ten vertices, which is simple, and forests are separated
by all of these invariants. Running `n = 2..9`:

| n | multigraphs with m ≤ 9 | collisions |
|---:|---:|---|
| 2–5 | 9 / 43 / 233 / 722 | none |
| 6 | 1 462 | one, at m = 9 |
| 7 | 1 738 | one, at m = 9 |
| 8 | 1 275 | one, at m = 9 |
| 9 | 526 | one, at m = 9 |

So **no counterexample has fewer than nine edges**, and nine is attained. The
six-vertex witness has maximum multiplicity two, so it also sits in the `(6,2)`
row above — it is the single `XB`-class in the `m = 9` bucket of the `(6,4)`
run. The five-vertex pair of the paper is smallest in *vertices*, not in edges.

Cost is dominated by computing `Aut(G)` over the `n!` permutations, so `n = 9`
takes a couple of minutes. `n = 10` is not attempted and is not needed.

## How it works

Connected-block polynomials `C_G(B;q)` are computed for all `2^n` vertex
subsets by the subset recursion (`O(3^n)` operations); each invariant is then a
set-partition sum over blocks. Expanding `C_G(B;q)` in the basis
`q^{s-1}(1+q)^r` is avoided by

```
sum_r a_r theta^r  =  C_G(B; theta-1) / (theta-1)^{s-1}
```

so a substitution `u_{s,r} -> sum_k beta_{k,s} theta_k^r` needs only `K`
numerical evaluations, and all three invariants share them. Arithmetic is in
`F_p`, `p = 2^61 - 1`, at random points. Fingerprint equality is a *necessary*
condition for equality of the invariants, so no collision can be missed; false
positives are removed afterwards by exact recomputation over **Z**.

Before any `O(3^n)` work, a staged sieve discards graphs using invariants that
`L_G` determines but that are much cheaper — `m`, the degree sequence and the
induced edge-count profile (`~2^n`); then spanning trees and triangles
(`~n^3`); then `X_G = XB_G(-1)` (`2*3^n`). On the complete `n = 9` list these
retain 9.2%, 0.18% and 0.058% respectively.

Every filter is determined by `L_G` **and** by `XB_G`, so a discarded graph
cannot belong to a collision class of either. Do not add filters known only for
`XB` — they could hide exactly the counterexample being searched for.

## Reproducibility notes

* Rerun with `-s SEED` for a second, independent set of substitution points;
  the classes must be identical.
* `--no-profile` drops the profile from the sieve: slower, but an independent
  path to the same answer.
* The `n <= 7` fingerprints were validated against `gp.py`, which implements
  the invariants directly from their defining expansions.
* For `n <= 8` the XB-classes reproduce those found by Markström for the
  U-polynomial by an unrelated method.

## Citing

If you use this code, please cite both the paper and the archived release:

> A. Kirillov, G. Nenashev, B. Shapiro and A. Vaintrob,
> *loopy-collisions: exhaustive collision search for the loopy polynomial*,
> Zenodo, 2026. doi:10.5281/zenodo.21929427

That is the *concept* DOI: it always resolves to the newest version.  To pin a
specific one, use its own DOI instead (v1.0.0 is
`10.5281/zenodo.21929428`).

`CITATION.cff` carries the same information in machine-readable form.

## License

MIT — see `LICENSE`.
