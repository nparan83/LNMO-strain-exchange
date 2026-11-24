#!/usr/bin/env python3
# pilot_dE.py -- LNMO_Pilots / LNMO_VC_Relax cikti taramasi
# paper-1'in dE(AFM-FM), dE(FiM-FM) degerlerinin kaynagini belirler
import re, glob, os
RY = 13605.693122994          # meV
FM40 = -6557.44511478         # Ry, 40 atomlu cfg0 (durum sec.7a)
FM20 = FM40/2.0

def parse(fn):
    t = open(fn, errors="ignore").read()
    e = re.findall(r"^!\s+total energy\s*=\s*(-?[\d.]+)", t, re.M)
    if not e: e = re.findall(r"^!!?\s*(-?\d+\.\d+)", t, re.M)
    nat = re.search(r"number of atoms/cell\s*=\s*(\d+)", t)
    mag = re.findall(r"total magnetization\s*=\s*(-?[\d.]+)", t)
    vol = re.findall(r"unit-cell volume\s*=\s*([\d.]+)", t)
    fin = "Final scf calculation" in t or "Final energy" in t
    nst = len(re.findall(r"number of bfgs steps", t)) or len(e)
    return (float(e[-1]) if e else None,
            int(nat.group(1)) if nat else None,
            float(mag[-1]) if mag else None,
            float(vol[-1])*0.14818471 if vol else None,   # bohr^3 -> A^3
            fin, nst, len(e))

files = sorted(set(glob.glob("LNMO_Pilots/*/*.out") + glob.glob("LNMO_VC_Relax/*.out")))
print("%-42s %5s %14s %7s %9s %6s %5s %12s" % ("dosya","nat","E(Ry)","Mtot","V(A^3)","final","nscf","dE(meV/f.u.)"))
for fn in files:
    E, nat, mag, vol, fin, nst, ne = parse(fn)
    if E is None: continue
    if nat in (20, None):  ref, nfu = FM20, 2.0
    elif nat == 40:        ref, nfu = FM40, 4.0
    else:                  ref, nfu = FM20*nat/20.0, nat/10.0
    dE = (E - ref)*RY/nfu
    print("%-42s %5s %14.8f %7s %9s %6s %5d %12.2f" % (
        fn[-42:], nat, E, ("%.2f"%mag) if mag is not None else "-",
        ("%.2f"%vol) if vol else "-", "E" if fin else "-", ne, dE))
print()
print("referans E_FM(20 atom) = %.8f Ry  (= cfg0/2, durum sec.7a)" % FM20)
print("paper-1 hedefleri: AFM 133.83 | FiM 192.2 (eps=0), 194.9 (-5%), 156.8 (+5%)")
print("paper-2 statik (FM geometrisinde): C-tipi 146.47 | G-tipi 215.96 (eps=0)")
