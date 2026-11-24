#!/usr/bin/env python3
# afm_id2.py -- QE girdisindeki spin desenini bag-bag A/C/G olarak teshis et
import re, sys, glob
import numpy as np
from q3lib import read_qe, neighbours, raw_labels, species


for fn in sorted(sys.argv[1:] or glob.glob("LNMO_Pilots/AFM_000/*.in")):
    sp, sm = species(fn); labs = raw_labels(fn)
    cell, els, X = read_qe(fn)
    if len(labs) != len(els): sys.exit("HATA: etiket sayisi (%d) atom sayisi (%d) ile uyusmuyor" % (len(labs), len(els)))
    sig = np.array([np.sign(sm.get(l, 0.0)) for l in labs])
    mag = [(l, sm.get(l, 0.0)) for l in sp]
    nz = [i for i in range(len(els)) if els[i] in ("Ni", "Mn") and sig[i] == 0]
    if nz: sys.exit("HATA: %d manyetik atomun isareti 0 -> etiket eslesmesi bozuk: %s" %
                    (len(nz), sorted({labs[i] for i in nz})))
    iNi = [i for i, s in enumerate(els) if s == "Ni"]
    iMn = [i for i, s in enumerate(els) if s == "Mn"]
    iO  = [i for i, s in enumerate(els) if s == "O"]
    pp = pa = qp = qa = 0
    for i in iNi:
        nb = neighbours(cell, X, i, iO, 6)
        nb.sort(key=lambda t: abs(t[1][2]), reverse=True)
        for k, (d, v, jo) in enumerate(nb):
            mn = neighbours(cell, np.vstack([X, X[i] + v]), len(X), iMn, 1)[0]
            same = sig[i] * sig[mn[2]] > 0
            if k < 2: qp += same; qa += not same
            else:     pp += same; pa += not same
    Mtot = sum(3.0 * sig[i] for i in iMn) + sum(2.0 * sig[i] for i in iNi)
    print("%s" % fn)
    print("  starting_magnetization:", ", ".join("%s=%+.2f" % m for m in mag))
    print("  duzlem-ici : %2d paralel / %2d antiparalel" % (pp, pa))
    print("  duzlem-disi: %2d paralel / %2d antiparalel" % (qp, qa))
    print("  nominal Mtot = %+.1f uB/hucre" % Mtot)
    tip = ("G-TIPI (FiM, hepsi AF)" if pa and qa else
           "C-TIPI (duzlem-ici AF)" if pa else
           "A-TIPI (duzlem-disi AF)" if qa else "FM")
    print("  -> %s" % tip)
