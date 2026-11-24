#!/usr/bin/env python3
# mktab.py -- Appendix tablolarini .dat dosyalarindan uretir (markdown + LaTeX).
# Hicbir sayi elle yazilmaz; her sutun bir dosya yoluna izlenir.
import os, sys, glob, math

D = sys.argv[1] if len(sys.argv) > 1 else "."

def load(name, need=None):
    p = os.path.join(D, name)
    if not os.path.exists(p):
        print("UYARI: bulunamadi ->", p); return None
    rows = []
    for l in open(p):
        s = l.strip()
        if not s or s.startswith("#"): continue
        t = s.split()
        v = []
        for x in t:
            try: v.append(float(x))
            except ValueError: v.append(None)
        if v and v[0] is not None:
            if need is None or len(v) >= need: rows.append(v)
    rows.sort(key=lambda r: r[0])
    return rows

def pct(rows):
    if not rows: return rows
    m = max(abs(r[0]) for r in rows)
    if m <= 1.0:
        for r in rows: r[0] *= 100.0
    return rows

def at0(rows, j):
    for r in rows:
        if abs(r[0]) < 1e-6: return r[j] if j < len(r) else None
    return None

def findcol(rows, target, tol, skip=1):
    z = None
    for r in rows:
        if abs(r[0]) < 1e-6: z = r
    if z is None: return None
    best, bj = None, None
    for j in range(skip, len(z)):
        if z[j] is None: continue
        d = abs(z[j] - target)
        if best is None or d < best: best, bj = d, j
    if best is None or best > tol: return None
    return bj

def se(e):
    return "0" if abs(e) < 1e-9 else "%+.0f" % e

def f(x, n):
    return "nan" if (x is None or (isinstance(x, float) and math.isnan(x))) else ("%.*f" % (n, x))

geo  = pct(load("q3_geo.dat", 11) or [])
dcD  = pct(load("dc_D.dat", 12) or [])
erel = pct(load("Erel.dat", 16) or [])
usc  = load("U_scan.dat", 8) or []
jstr = pct(load("J_strain.dat", 3) or [])
tcs  = pct(load("Tc_strain.dat", 2) or [])

print("okunan satir sayilari: q3_geo=%d dc_D=%d Erel=%d U_scan=%d J_strain=%d Tc_strain=%d"
      % (len(geo), len(dcD), len(erel), len(usc), len(jstr), len(tcs)))

REF = [("q3_geo th_par(0)",   at0(geo, 9),  155.780, 0.01),
       ("q3_geo th_perp(0)",  at0(geo, 10), 155.156, 0.01),
       ("q3_geo Q3_Ni(0)",    at0(geo, 3),  -0.003469, 0.002),
       ("dc_D c_perp(0)",     at0(dcD, 3),  7.845544, 1e-4),
       ("dc_D beta(0)",       at0(dcD, 5),  125.192, 0.01),
       ("dc_D a_pc(0)",       at0(dcD, 6),  3.93476, 1e-3),
       ("dc_D J_par(0)",      at0(dcD, 9),  12.206, 0.01),
       ("dc_D J_perp(0)",     at0(dcD, 10), 11.582, 0.01),
       ("dc_D D_par/D_perp(0)", at0(dcD,11), 1.0603, 0.002),
       ("Erel V_FM(0)",       at0(erel, 2), 242.92, 0.05),
       ("Erel dE_AFM(0)",     at0(erel, 5), 133.83, 0.05),
       ("Erel Erel_C(0)",     at0(erel, 7), 12.64,  0.05)]
bad = 0
for nm, got, want, tol in REF:
    ok = (got is not None) and abs(got - want) <= tol
    if not ok: bad += 1
    print("  %-24s got=%-12s ref=%-10s %s" % (nm, f(got, 5), want, "OK" if ok else "<-- SAPMA"))
print("dogrulama: %d/%d sapma" % (bad, len(REF)))

jp = findcol(jstr, 12.206, 0.02) if jstr else None
jq = findcol(jstr, 11.582, 0.02) if jstr else None
tc = findcol(tcs, 299.9, 1.5) if tcs else None
print("J_strain: J_par sutunu=%s  J_perp sutunu=%s | Tc_strain: T_C sutunu=%s" % (jp, jq, tc))

L1 = []
L1.append("**Table 1.** Relaxed geometry of the strain series. Sources: `q3_geo.dat` "
          "(distances, angles, Q3), `dc_D.dat` (lattice), `Erel.dat` (volume).\n")
L1.append("| eps (%) | d(Ni-O)_eq (A) | d(Ni-O)_ap (A) | th_par (deg) | th_perp (deg) "
          "| Q3(Ni) (A) | Q3(Mn) (A) | c_perp (A) | beta (deg) | a_pc (A) | V (A^3) |")
L1.append("|" + "---|" * 11)
for i, r in enumerate(geo):
    e = r[0]
    dd = next((q for q in dcD if abs(q[0] - e) < 1e-6), None)
    ee = next((q for q in erel if abs(q[0] - e) < 1e-6), None)
    L1.append("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
        se(e), f(r[1], 4), f(r[2], 4), f(r[9], 3), f(r[10], 3), f(r[3], 5), f(r[7], 5),
        f(dd[3] if dd else None, 5), f(dd[5] if dd else None, 3),
        f(dd[6] if dd else None, 5), f(ee[2] if ee else None, 2)))

L2 = []
L2.append("\n**Table 2.** Exchange constants. Source: `dc_D.dat` (J_par, J_perp, D ratio)"
          "%s.\n" % ("; `J_strain.dat` (cross-check)" if jp is not None else ""))
L2.append("| eps (%) | J_par (meV) | J_perp (meV) | J_par/J_perp | J_par-J_perp (meV) "
          "| <J1> (meV) | D_par/D_perp |")
L2.append("|" + "---|" * 7)
for r in dcD:
    a, b = r[9], r[10]
    L2.append("| %s | %s | %s | %s | %s | %s | %s |" % (
        se(r[0]), f(a, 3), f(b, 3), f(a / b, 4), f(a - b, 3), f((2 * a + b) / 3.0, 3), f(r[11], 4)))
if jp is not None and jq is not None:
    dm = max(abs(r[jp] - next(q[9] for q in dcD if abs(q[0] - r[0]) < 1e-6)) for r in jstr
             if any(abs(q[0] - r[0]) < 1e-6 for q in dcD))
    L2.append("\n*(J_strain.dat vs dc_D.dat J_par: max fark %.4f meV)*" % dm)

L3 = []
if tc is not None:
    L3.append("\n**Table 3.** Monte-Carlo Curie temperatures. Source: `Tc_strain.dat` "
              "(from `mc_eps*.dat`).\n")
    L3.append("| eps (%) | T_C (K) | T_C/T_C(0) |")
    L3.append("|---|---|---|")
    t0 = at0(tcs, tc)
    for r in tcs:
        L3.append("| %s | %s | %s |" % (se(r[0]), f(r[tc], 1), f(r[tc] / t0, 3)))
else:
    L3.append("\n*(Table 3: Tc_strain.dat okunamadi veya T_C sutunu bulunamadi.)*")

L4 = []
L4.append("\n**Table 4.** Hubbard-U scan. Source: `U_scan.dat`.\n")
L4.append("| eps (%) | U(Mn) | U(Ni) | J_par | J_perp | J_par/J_perp | J_par-J_perp | <J1> | max resid |")
L4.append("|" + "---|" * 9)
for r in sorted(usc, key=lambda x: (x[0], x[2])):
    e = r[0] * 100.0 if abs(r[0]) <= 1.0 and r[0] != 0 else r[0]
    L4.append("| %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
        se(e), f(r[1], 1), f(r[2], 1), f(r[3], 4), f(r[4], 4), f(r[5], 4),
        f(r[3] - r[4], 4), f(r[6], 4), f(r[7], 3)))

out = "\n".join(L1 + L2 + L3 + L4) + "\n"
open(os.path.join(D, "tables.md"), "w").write(out)

def tex(lines):
    body = [l for l in lines if l.startswith("|") and not l.startswith("|---")]
    if not body: return ""
    hdr = [c.strip() for c in body[0].strip("|").split("|")]
    n = len(hdr)
    s = ["\\begin{tabular}{" + "r" * n + "}", "\\hline",
         " & ".join(h.replace("_", "\\_").replace("%", "\\%").replace("<", "$\\langle$")
                     .replace(">", "$\\rangle$") for h in hdr) + " \\\\", "\\hline"]
    for l in body[1:]:
        s.append(" & ".join(c.strip() for c in l.strip("|").split("|")) + " \\\\")
    s += ["\\hline", "\\end{tabular}", ""]
    return "\n".join(s)

open(os.path.join(D, "tables.tex"), "w").write(
    tex(L1) + "\n" + tex(L2) + "\n" + tex(L3) + "\n" + tex(L4))
print("yazildi: tables.md (%d satir) + tables.tex" % len(out.splitlines()))
