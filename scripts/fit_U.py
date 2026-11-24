#!/usr/bin/env python3
# fit_U.py -- U taramasi sonucu: J_par(U), J_perp(U), anizotropi dayanikliligi
import os, re, sys, glob
import numpy as np
RY = 13605.693122994
USETS = [("base", 4.0, 5.1), ("m10", 3.0, 4.1), ("p10", 5.0, 6.1)]
STRAINS = [-5, 0, 5]; CFGS = ["FM", "A", "C", "G"]
# 40 atomlu hucre, sec.5b: sutunlar [1, J_par, J_perp]
M = np.array([[1.,-24.,-12.],[1.,-24.,+12.],[1.,+24.,-12.],[1.,+24.,+12.]])

def energy(fn):
    if not os.path.exists(fn): return None
    t = open(fn, errors="ignore").read()
    e = re.findall(r"^!\s+total energy\s*=\s*(-?[\d.]+)", t, re.M)
    ok = "JOB DONE" in t or "convergence has been achieved" in t
    return (float(e[-1]), ok) if e else None

def fname(tag, e, c):
    if tag == "base":
        for p in ("e%+03d_%s.out" % (e, c), "e%02d_%s.out" % (e, c),
                  {"FM":"eps0_cfg0.out","A":"eps0_cfgA.out","C":"eps0_cfgC.out","G":"eps0_cfg1.out"}[c]):
            if os.path.exists(p): return p
        return "e%+03d_%s.out" % (e, c)
    return "U%s_e%+03d_%s.out" % (tag, e, c)

rows, miss = [], []
for tag, uMn, uNi in USETS:
    for e in STRAINS:
        E, bad = [], False
        for c in CFGS:
            r = energy(fname(tag, e, c))
            if r is None: miss.append(fname(tag, e, c)); bad = True; break
            if not r[1]: print("UYARI: yakinsama supheli ->", fname(tag, e, c))
            E.append(r[0])
        if bad: continue
        y = (np.array(E) - E[0]) * RY
        p, *_ = np.linalg.lstsq(M, y, rcond=None)
        res = y - M @ p; Jp, Jq = p[1], p[2]
        rows.append([e, tag, uMn, uNi, Jp, Jq, Jp/Jq, (4*Jp+2*Jq)/6, np.abs(res).max()])

if miss:
    print("EKSIK cikti (%d):" % len(miss)); [print("  ", m) for m in miss[:8]]
if not rows: sys.exit("hicbir set tamamlanmamis")
print("%5s %6s %5s %5s %9s %9s %9s %9s %9s" % ("eps","Uset","U_Mn","U_Ni","J_par","J_perp","J||/Jp","<J1>","|artik|"))
for r in rows: print("%+5.0f %6s %5.1f %5.1f %9.4f %9.4f %9.4f %9.4f %9.3f" % tuple(r))
with open("U_scan.dat", "w") as f:
    f.write("# eps  U_Mn  U_Ni  J_par  J_perp  ratio  <J1>  max|artik|\n")
    for r in rows: f.write("%+5.0f %5.1f %5.1f %10.5f %10.5f %10.5f %10.5f %9.3f\n" % (r[0], *r[2:]))
print()
for e in STRAINS:
    S = [r for r in rows if r[0] == e]
    if len(S) < 2: continue
    u = np.array([s[3] for s in S])
    for nm, k in (("J_par",4), ("J_perp",5), ("J||/Jperp",6)):
        v = np.array([s[k] for s in S]); c = np.polyfit(u, v, 1)
        print("eps=%+d  d(%s)/dU = %+.4f / eV   (%.4f -> %.4f, %%%+.1f)" % (
              e, nm, c[0], v[np.argmin(u)], v[np.argmax(u)], 100*(v[np.argmax(u)]/v[np.argmin(u)]-1)))
    print()
if len({r[0] for r in rows}) >= 2:
    print("--- ANA SORU: anizotropinin strain tepkisi U'ya dayanikli mi? ---")
    for tag in [t for t, _, _ in USETS]:
        a = [r for r in rows if r[1] == tag and r[0] == 0]
        b = [r for r in rows if r[1] == tag and r[0] == 5]
        if a and b:
            print("  U=%-5s : J||/Jperp  %.4f (eps=0) -> %.4f (eps=+5)   degisim %+.4f" % (
                  tag, a[0][6], b[0][6], b[0][6]-a[0][6]))
    print("  (taban: 1.0539 -> 1.5497, degisim +0.4958)")
