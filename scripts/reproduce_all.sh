#!/bin/bash
# reproduce_all.sh -- regenerate every derived data file in data/ from the Quantum ESPRESSO outputs and compare.
# Usage:  bash scripts/reproduce_all.sh [path to a clone of LNMO-strain-parity]
# Without the argument the geometry files (q3_geo.dat, dc_D.dat, Erel.dat) are taken from data/; with it they are
# rebuilt from the relaxed structures and vc-relax outputs of the companion data set (doi:10.5281/zenodo.24353930).
# Work directory: ./work_reproduce (or $WORK). Requires Python >= 3.8 with NumPy. Run time: < 1 min.
set -eu
REPO=$(cd "$(dirname "$0")/.." && pwd); P1=${1:-}; [ -n "$P1" ] && P1=$(cd "$P1" && pwd)
W=${WORK:-$PWD/work_reproduce}; rm -rf "$W"; mkdir -p "$W"; cd "$W"
for d in strain_2x2 U_scan; do for f in "$REPO"/outputs/$d/*.out.gz; do gunzip -c "$f" > "$(basename "${f%.gz}")"; done; done
cp "$REPO"/scripts/*.py .; cp "$REPO"/data/mc/mc_eps*.dat .
python3 fit_all.py > fit_all.log
if [ -n "$P1" ]; then
  for d in "$P1"/inputs/strain_*; do e=$(basename "$d"); for c in FM AFM FiM; do [ -d "$d/$c" ] || continue
    t=LNMO_Strain_Production/$e/$c; mkdir -p "$t"; cp "$d/$c/scf_$c.in" "$t/"
    for f in "$P1/outputs/$e/$c"/scf_FM.out.gz "$P1/outputs/$e/$c"/scf_final_*.out.gz "$P1/outputs/$e/$c"/vcrelax_*.out.gz; do
      [ -f "$f" ] && gunzip -c "$f" > "$t/$(basename "${f%.gz}")"; done; done; done
  python3 q3_geo.py LNMO_Strain_Production > q3_geo.log
  python3 q3_align.py > q3_align.log; python3 getc.py > getc.log; python3 mkerel.py > mkerel.log
else
  cp "$REPO"/data/q3_geo.dat "$REPO"/data/dc_D.dat "$REPO"/data/Erel.dat .; python3 q3_align.py > q3_align.log
fi
python3 tc.py > tc.log; python3 fit_U.py > fit_U.log; python3 mktab.py . > mktab.log
python3 - "$REPO/data" << 'PY'
import sys, os, numpy as np
R = sys.argv[1]; bad = 0
def rows(f, drop_last=False):
    out = []
    for l in open(f):
        if l.startswith('#') or not l.strip(): continue
        s = l.split(); out.append([float(x) for x in (s[:-1] if drop_last else s)])
    a = np.array(out); return a[np.argsort(a[:, 0])]
for f in ('J_strain.dat', 'q3_geo.dat', 'q3_align.dat', 'dc_D.dat', 'Erel.dat', 'Tc_strain.dat', 'U_scan.dat'):
    a, b = rows(f, f == 'q3_geo.dat'), rows(os.path.join(R, f), f == 'q3_geo.dat')
    d = np.abs(a - b); d[np.isnan(d)] = 0; ok = a.shape == b.shape and d.max() < 1e-6; bad += not ok
    print('%-14s %-8s max|diff| = %.1e' % (f, 'OK' if ok else 'DIFFERS', d.max()))
for f in ('tables.md', 'tables.tex'):
    ok = open(f).read() == open(os.path.join(R, f)).read(); bad += not ok
    print('%-14s %s' % (f, 'OK (identical)' if ok else 'DIFFERS'))
print('result: %d file(s) differ' % bad)
PY
