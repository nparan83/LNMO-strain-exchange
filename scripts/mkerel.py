#!/usr/bin/env python3
# mkerel.py -- manyetostruktürel gevsetme enerjisi E_rel(eps) + FM elastik egrisi
# Girdi : LNMO_Strain_Production/strain_*/{FM,AFM,FiM}/*.out  +  q3_align.dat
# Cikti : Erel.dat  +  ekrana ~24 satir.  SLURM GEREKMEZ.
import os, re, sys, glob
import numpy as np

ROOT = sys.argv[1] if len(sys.argv) > 1 else "LNMO_Strain_Production"
ALG  = sys.argv[2] if len(sys.argv) > 2 else "q3_align.dat"
RY   = 13605.693122994          # meV
BOHR3 = 0.14818471              # (a.u.)^3 -> A^3

def eps_of(p):
    m = re.search(r"strain_([+-]?[pm]?\d+(?:[._]\d+)?)", p)
    if not m: return None
    return float(m.group(1).replace("_", ".").replace("p", "+").replace("m", "-"))

def parse(fn):
    t = open(fn, errors="ignore").read()
    e = re.findall(r"^!\s+total energy\s*=\s*(-?[\d.]+)", t, re.M)
    if not e: return None
    v = re.findall(r"unit-cell volume\s*=\s*([\d.]+)", t)
    n = re.search(r"number of atoms/cell\s*=\s*(\d+)", t)
    return dict(E=float(e[-1]), V=float(v[-1])*BOHR3 if v else np.nan,
                nat=int(n.group(1)) if n else 0,
                relax=("bfgs" in t or "Final scf" in t), f=os.path.basename(fn))

# oncelik: son-scf > vcrelax > diger.  Ayni dizindeki TUM adaylar karsilastirilir.
def pick(d):
    cand = [parse(f) for f in sorted(glob.glob(os.path.join(d, "*.out")))]
    cand = [c for c in cand if c]
    if not cand: return None, []
    def rank(c):
        n = c["f"]
        return (0 if n.startswith("scf_final") else 1 if n.startswith("scf_") else
                2 if n.startswith("vcrelax") else 3, n)
    cand.sort(key=rank)
    return cand[0], cand

rows, notes = [], []
for d in sorted(glob.glob(os.path.join(ROOT, "strain_*")), key=lambda p: eps_of(p)):
    e = eps_of(d)
    r = {"eps": e*100 if abs(e) <= 1 else e}
    for cfg in ("FM", "AFM", "FiM"):
        best, cand = pick(os.path.join(d, cfg))
        r[cfg] = best
        if len(cand) > 1:
            dE = (max(c["E"] for c in cand) - min(c["E"] for c in cand)) * RY / 2
            if dE > 0.5:
                notes.append("eps=%+.0f %s: adaylar %.2f meV/f.u. farkli -> %s" % (
                    r["eps"], cfg, dE, ", ".join("%s(%.6f)" % (c["f"], c["E"]) for c in cand)))
    rows.append(r)
if not rows: sys.exit("HATA: dosya bulunamadi -> " + ROOT)

T = np.loadtxt(ALG, usecols=range(8)); T = T[np.argsort(T[:, 0])]
if np.abs(T[:, 0]).max() <= 1.0: T[:, 0] *= 100.0
JP = {int(round(a)): (b, c) for a, b, c in zip(T[:, 0], T[:, 4], T[:, 5])}

out = []
for r in rows:
    e = r["eps"]; fm = r["FM"]
    if fm is None: sys.exit("HATA: eps=%+.0f icin FM ciktisi yok" % e)
    jp, jq = JP[int(round(e))]
    Cst = 12.0*jp                      # dE(C-FM) statik, meV/f.u.
    Gst = 1.5*(8.0*jp + 4.0*jq)        # dE(G-FM) statik, meV/f.u.
    def dE(c): return (c["E"] - fm["E"])*RY/2.0 if c else np.nan
    af, fi = r["AFM"], r["FiM"]
    out.append([e, fm["E"], fm["V"],
                af["E"] if af else np.nan, af["V"] if af else np.nan, dE(af), Cst, Cst-dE(af),
                (af["V"]-fm["V"]) if af else np.nan,
                fi["E"] if fi else np.nan, fi["V"] if fi else np.nan, dE(fi), Gst, Gst-dE(fi),
                (fi["V"]-fm["V"]) if fi else np.nan])
O = np.array(out, float)
i0 = int(np.argmin(np.abs(O[:, 0])))
EL = (O[:, 1] - O[i0, 1])*RY/2.0                      # FM elastik egrisi, meV/f.u.
O = np.column_stack([O, EL])

np.savetxt("Erel.dat", O, fmt="%+13.6f",
  header=("eps E_FM(Ry) V_FM E_AFM V_AFM dE_AFM C_stat Erel_C dV_AFM "
          "E_FiM V_FiM dE_FiM G_stat Erel_G dV_FiM E_elastic\n"
          "meV/f.u. birimleri; V: A^3 (20 atomlu hucre); nan = o noktada konfig yok\n"
          "kaynak: %s + %s" % (ROOT, ALG)))

print("Erel.dat yazildi: %d nokta" % len(O))
print("secilen dosyalar (eps=0): FM=%s  AFM=%s  FiM=%s" % (
      rows[i0]["FM"]["f"], rows[i0]["AFM"]["f"] if rows[i0]["AFM"] else "-",
      rows[i0]["FiM"]["f"] if rows[i0]["FiM"] else "-"))
for n in notes[:4]: print("  ** " + n)
if not notes: print("  aday dosyalar arasi enerji farki < 0.5 meV/f.u. -> secim guvenli")
print()
print("%5s %9s %9s %8s %8s | %9s %9s %8s %8s" % (
      "eps","dE_AFM","C_stat","Erel_C","dV_AFM","dE_FiM","G_stat","Erel_G","dV_FiM"))
for r in O:
    print("%+5.0f %9.2f %9.2f %8.2f %8.3f | %9s %9.2f %8s %8s" % (
        r[0], r[5], r[6], r[7], r[8],
        "%.2f" % r[11] if np.isfinite(r[11]) else "  -  ", r[12],
        "%.2f" % r[13] if np.isfinite(r[13]) else "  -  ",
        "%+.3f" % r[14] if np.isfinite(r[14]) else "  -  "))
m = np.isfinite(O[:, 7])
print("\nE_rel(C) = %.2f +- %.2f meV/f.u.  (n=%d, tepe eps=%+.0f)" % (
      O[m,7].mean(), O[m,7].std(ddof=1), m.sum(), O[np.argmax(np.where(m,O[:,7],-9e9)),0]))
g = np.isfinite(O[:, 13])
print("E_rel(G) = %s meV/f.u.  (n=%d)" % (" / ".join("%.2f" % v for v in O[g,13]), g.sum()))
print("FM elastik: -5 -> %.1f ; +5 -> %.1f meV/f.u.  asimetri %.1f" % (EL[0], EL[-1], EL[0]-EL[-1]))
c = np.polyfit(O[:,0], EL, 2)
print("            parabol %.3f eps^2 %+.3f eps %+.3f  min @ eps=%+.3f%%" % (*c, -c[1]/(2*c[0])))
REFC = {-5:11.52,-4:11.71,-3:12.29,-2:12.54,-1:12.52,0:12.65,1:12.13,2:12.18,3:11.97,4:11.94,5:11.79}
REFE = {-5:552.7, 0:0.0, 5:428.0}
bad = [(int(r[0]), round(r[7],2), REFC[int(r[0])]) for r in O
       if int(r[0]) in REFC and abs(r[7]-REFC[int(r[0])]) > 0.05]
bad += [("elast", round(EL[i],1), REFE[int(O[i,0])]) for i in range(len(O))
        if int(O[i,0]) in REFE and abs(EL[i]-REFE[int(O[i,0])]) > 0.5]
print("sec.7k karsilastirmasi:", "DOGRULANDI" if not bad else "*** SAPMA %s" % bad)
