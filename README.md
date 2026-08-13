# loopy-collisions

Code accompanying **A. Kirillov, G. Nenashev, B. Shapiro and A. Vaintrob,
*The Loopy Polynomial: from Tutte's Universal V-Function to Bizonotopal
Geometry*** (Section 4.3, "The computation").

It searches exhaustively for pairs of non-isomorphic connected simple graphs
with equal **loopy polynomial** `L_G`, equal **refined loopy polynomial**
`Lhat_G`, or equal **Tutte symmetric function** `XB_G`, and checks that the
three induce the same partition. Verified through **n = 11** vertices
(1 006 700 565 connected graphs).

## Quick start

```sh
make            # builds loopy_collisions, comp_check, gen_graphs
make check      # reproduces the published n <= 9 answers, a few seconds
```

`make check` should print, for `n = 9`:

```
SUMMARY L n=9 graphs=152 classes=65 pairs=65 edges=[13,23]
VERDICT  L-classes == Lhat-classes : YES
VERDICT  Lhat-classes == XB-classes: YES
VERDICT  L-classes == XB-classes   : YES
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
> version 1.0.0, Zenodo, 2026. doi:10.5281/zenodo.XXXXXXX

`CITATION.cff` carries the same information in machine-readable form.

## License

MIT — see `LICENSE`.
