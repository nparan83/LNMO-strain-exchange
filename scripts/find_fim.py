#!/usr/bin/env python3
# find_fim.py -- Mtot ~ 2 (FiM) imzasi tasiyan tum QE ciktilarini bul
import os, re
RY = 13605.693122994
FM20 = -6557.44511478/2.0
rows = []
for root, dirs, fs in os.walk("."):
    dirs[:] = [d for d in dirs if not (d.endswith(".save") or d == "tmp" or d == "out")]
    for fn in fs:
        if not fn.endswith(".out"): continue
        p = os.path.join(root, fn)
        try:
            if os.path.getsize(p) > 300*1024*1024: continue
            t = open(p, errors="ignore").read()
        except Exception: continue
        e = re.findall(r"^!\s+total energy\s*=\s*(-?[\d.]+)", t, re.M)
        if not e: continue
        m = re.findall(r"total magnetization\s*=\s*(-?[\d.]+)", t)
        n = re.search(r"number of atoms/cell\s*=\s*(\d+)", t)
        v = re.findall(r"unit-cell volume\s*=\s*([\d.]+)", t)
        rows.append((p, int(n.group(1)) if n else 0, float(e[-1]),
                     float(m[-1]) if m else None,
                     float(v[-1])*0.14818471 if v else None,
                     "vcrelax" if "bfgs" in t or "Final scf" in t else "scf", len(e)))
print("toplam taranan .out (enerjili):", len(rows))
h = {}
for r in rows:
    k = "None" if r[3] is None else "%.1f" % r[3]
    h[k] = h.get(k, 0) + 1
print("Mtot histogrami:", " | ".join("%s:%d" % kv for kv in sorted(h.items(), key=lambda x: x[1], reverse=True)[:10]))
print()
print("--- Mtot ~ 2 (FiM adaylari) ---")
print("%-46s %4s %15s %6s %8s %8s %11s" % ("dosya","nat","E(Ry)","Mtot","V(A^3)","tip","dE(meV/fu)"))
hit = [r for r in rows if r[3] is not None and 1.5 < abs(r[3]) < 2.5]
for p, nat, E, m, v, typ, ne in sorted(hit):
    nfu = max(nat/10.0, 1e-9); ref = FM20*nat/20.0
    print("%-46s %4d %15.8f %6.2f %8s %8s %11.2f" % (
        p[-46:], nat, E, m, ("%.2f"%v) if v else "-", typ, (E-ref)*RY/nfu))
if not hit: print("  (yok)")
print()
print("--- eski pilot serisi lmno_str_* ---")
for p, nat, E, m, v, typ, ne in sorted(rows):
    if "lmno_str" in p:
        nfu = max(nat/10.0, 1e-9); ref = FM20*nat/20.0
        print("%-46s %4d %15.8f %6s %8s %11.2f" % (p[-46:], nat, E,
              ("%.2f"%m) if m is not None else "-", ("%.2f"%v) if v else "-", (E-ref)*RY/nfu))
print()
print("paper-1 FiM hedefleri: 192.2 (eps=0) | 194.9 (-5%) | 156.8 (+5%)")
