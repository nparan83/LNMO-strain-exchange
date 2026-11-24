#!/usr/bin/env python3
# getc.py -- 11 strain noktasinda gercek c_perp, beta, a_pc ve D_par/D_perp
# Girdi : LNMO_Strain_Production/strain_*/FM/scf_FM.in  +  q3_align.dat
# Cikti : dc_D.dat  +  ekrana 14 satir.  SLURM GEREKMEZ.
import os, re, sys, glob
import numpy as np
sys.path.insert(0, ".")
from q3lib import read_qe

ROOT = sys.argv[1] if len(sys.argv) > 1 else "LNMO_Strain_Production"
ALG  = sys.argv[2] if len(sys.argv) > 2 else "q3_align.dat"

def eps_of(p):
    m = re.search(r"strain_([+-]?[pm]?\d+(?:[._]\d+)?)", p)
    if not m: return None
    return float(m.group(1).replace("_", ".").replace("p", "+").replace("m", "-"))

rows = []
for fn in sorted(glob.glob(os.path.join(ROOT, "strain_*", "FM", "scf_FM.in"))):
    e = eps_of(fn)
    cell, els, X = read_qe(fn)
    a1, a2, a3 = cell
    a, b = np.linalg.norm(a1), np.linalg.norm(a2)
    cperp = abs(a3[2]); cmag = np.linalg.norm(a3)
    beta = np.degrees(np.arccos(a1 @ a3 / (a * cmag)))
    rows.append([e, a, b, cperp, cmag, beta])
R = np.array(rows, float)
if len(R) == 0: sys.exit("HATA: dosya bulunamadi -> " + ROOT)
if np.abs(R[:, 0]).max() <= 1.0: R[:, 0] *= 100.0
R = R[np.argsort(R[:, 0])]
eps, a, b, cperp, cmag, beta = R.T
a_pc = (a + b) / (2 * np.sqrt(2)); c_pc = cperp / 2.0

T = np.loadtxt(ALG, usecols=range(8))
T = T[np.argsort(T[:, 0])]
if np.abs(T[:, 0]).max() <= 1.0: T[:, 0] *= 100.0
if len(T) != len(R) or np.abs(T[:, 0] - eps).max() > 1e-6:
    sys.exit("HATA: q3_align.dat ile eps eslesmiyor")
Jp, Jq, Q3 = T[:, 4], T[:, 5], T[:, 1]
DD = (Jp / Jq) * (a_pc / c_pc) ** 2

np.savetxt("dc_D.dat", np.column_stack([eps, a, b, cperp, cmag, beta, a_pc, c_pc, Q3, Jp, Jq, DD]),
           fmt="%+10.5f",
           header="eps  a  b  c_perp  |a3|  beta(deg)  a_pc  c_pc  Q3_Ni  J_par  J_perp  D_par/D_perp\n"
                  "kaynak: LNMO_Strain_Production/strain_*/FM/scf_FM.in + q3_align.dat")

print("dc_D.dat yazildi:", len(R), "nokta")
print("%6s %10s %10s %10s %10s %9s %9s" % ("eps","a","b","c_perp","beta","a_pc","D||/Dperp"))
for i in range(len(R)):
    print("%+6.1f %10.6f %10.6f %10.6f %10.4f %9.5f %9.4f" % (eps[i],a[i],b[i],cperp[i],beta[i],a_pc[i],DD[i]))
ref = {-5: 8.208995, -1: 7.909111, 0: 7.845544, 1: 7.786885, 5: 7.586216}
bad = [(k, v, cperp[int(round(k)) + 5]) for k, v in ref.items() if abs(cperp[int(round(k)) + 5] - v) > 1e-4]
print("sec.7e c_perp kontrolu:", "5/5 DOGRULANDI" if not bad else "*** SAPMA: %s" % bad)
i0 = int(np.argmin(np.abs(eps)))
sgn = DD - 1.0; kk = np.where(np.diff(np.sign(sgn)))[0]
iso = (eps[kk[0]] + (0 - sgn[kk[0]]) / (sgn[kk[0]+1] - sgn[kk[0]]) * (eps[kk[0]+1] - eps[kk[0]])) if len(kk) else float("nan")
print("D||/Dperp: %.4f (eps=-5) -> %.4f (eps=0) -> %.4f (eps=+5)   toplam %.2f kat" % (DD[0],DD[i0],DD[-1],DD[-1]/DD[0]))
print("D||=Dperp noktasi: eps = %+.3f%%" % iso)
print("a_pc(eps=0) = %.4f A   (paper-1: 3.935)" % a_pc[i0])
