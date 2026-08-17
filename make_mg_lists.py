#!/usr/bin/env python3
"""make_mg_lists.py -- collect a multigraph run into one canonical file.

    python3 make_mg_lists.py res_6_4 -o data/multigraphs_n6_mu4.txt

`run_multi.sh` leaves its output scattered over one file per edge bucket, with
the classes of all four invariants interleaved.  This script gathers them,
checks the internal consistency of the run, and writes a single sorted file:

  * a header with the per-bucket counts and the totals;
  * every nontrivial L-class and every nontrivial XB-class, sorted by
    (invariant, m, size, members).

It also re-derives the four VERDICT lines from the raw CLASS lines rather than
trusting the ones the C program printed, and exits nonzero if any of the
must-hold ones fails:

    L-classes == Lhat-classes      (Problem: does L determine Lhat?)
    XB-classes == U-classes        (Sarmiento's equivalence)

The third, Lhat-classes == XB-classes, is the one expected to FAIL on
multigraphs -- that failure is the point of the run -- so it is reported but
not treated as an error.
"""
import argparse
import os
import re
import sys
from collections import defaultdict

CLASS = re.compile(r'^CLASS\s+(\S+)\s+size=(\d+)\s+m=(\d+)\s*:\s*(.*)$')
SUMM = re.compile(r'^SUMMARY\s+(\S+)\s+n=(\d+)\s+graphs=(\d+)\s+classes=(\d+)')
MEMB = re.compile(r'\[([^\]]*)\]')

KINDS = ('L', 'Lhat', 'XB', 'U')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('resdir')
    ap.add_argument('-o', '--out', required=True)
    args = ap.parse_args()

    files = sorted((int(re.match(r'm(\d+)\.txt$', f).group(1)), f)
                   for f in os.listdir(args.resdir)
                   if re.match(r'm(\d+)\.txt$', f))
    if not files:
        sys.exit('no m*.txt files in %s' % args.resdir)

    classes = {k: defaultdict(list) for k in KINDS}   # kind -> m -> [members]
    graphs = {}                                       # m -> multigraphs read
    declared = defaultdict(dict)                      # m -> kind -> classes
    n = None

    for m, f in files:
        for line in open(os.path.join(args.resdir, f)):
            mo = CLASS.match(line)
            if mo:
                kind, _, mm, rest = mo.groups()
                if kind not in classes:
                    sys.exit('unknown invariant %r in %s' % (kind, f))
                if int(mm) != m:
                    sys.exit('bucket %s contains m=%s' % (f, mm))
                members = tuple(sorted(' '.join(g.split())
                                       for g in MEMB.findall(rest)))
                classes[kind][m].append(members)
                continue
            mo = SUMM.match(line)
            if mo:
                kind, nn, g, c = mo.groups()
                n = int(nn) if n is None else n
                if int(nn) != n:
                    sys.exit('mixed vertex counts: %d and %s' % (n, nn))
                graphs[m] = int(g)
                declared[m][kind] = int(c)

    # ---- consistency of the run itself -------------------------------------
    for m in sorted(declared):
        for kind in KINDS:
            got = len(classes[kind].get(m, []))
            want = declared[m].get(kind)
            if want is None:
                sys.exit('bucket m=%d has no SUMMARY line for %s' % (m, kind))
            if got != want:
                sys.exit('bucket m=%d: %d %s-classes found, SUMMARY says %d'
                         % (m, got, kind, want))
        for kind in KINDS:
            if len(set(classes[kind].get(m, []))) != len(classes[kind].get(m, [])):
                sys.exit('bucket m=%d: duplicate %s-classes' % (m, kind))

    def partition(kind):
        return {m: set(classes[kind].get(m, [])) for m in sorted(declared)}

    def same(a, b):
        return partition(a) == partition(b)

    ok = True
    verdicts = []
    for a, b, must in (('L', 'Lhat', True), ('XB', 'U', True),
                       ('Lhat', 'XB', False)):
        eq = same(a, b)
        verdicts.append('# VERDICT  %-4s-classes == %-4s-classes : %s%s'
                        % (a, b, 'YES' if eq else 'NO',
                           '' if (eq or not must) else '   <-- MUST HOLD'))
        if must and not eq:
            ok = False
        if not must and not eq:
            verdicts[-1] += '   <-- counterexample, as expected'

    # ---- write ------------------------------------------------------------
    tot = sum(graphs.values())
    totc = {k: sum(len(v) for v in classes[k].values()) for k in KINDS}
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, 'w') as fh:
        fh.write('# loopless multigraphs on n = %d vertices\n' % n)
        fh.write('# %d multigraphs; %s\n'
                 % (tot, ', '.join('%d %s-classes' % (totc[k], k) for k in KINDS)))
        fh.write('#\n')
        for line in verdicts:
            fh.write(line + '\n')
        fh.write('#\n# m      multigraphs   L   Lhat     XB      U\n')
        for m in sorted(declared):
            fh.write('# %-5d %11d %5d %6d %6d %6d\n'
                     % (m, graphs[m], len(classes['L'].get(m, [])),
                        len(classes['Lhat'].get(m, [])),
                        len(classes['XB'].get(m, [])),
                        len(classes['U'].get(m, []))))
        fh.write('#\n# One class per line:  INVARIANT  m  [edge:mult ...] ...\n#\n')
        for kind in ('L', 'XB'):
            for m in sorted(classes[kind]):
                for members in sorted(classes[kind][m]):
                    fh.write('%-4s %3d %s\n'
                             % (kind, m, ' '.join('[%s]' % g for g in members)))

    print('%s: %d multigraphs, %s'
          % (args.out, tot, ', '.join('%d %s' % (totc[k], k) for k in KINDS)))
    for line in verdicts:
        print(line[2:])
    if not ok:
        sys.exit('a must-hold verdict failed')


if __name__ == '__main__':
    main()
