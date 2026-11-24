#!/usr/bin/env python3
# mkU.py -- U duyarlilik taramasi girdi ureteci
# Mevcut 2^2 konfig girdilerini alir, SADECE HUBBARD U ve prefix/outdir degistirir.
# Geometriye, k-orgusune, hicbir seye dokunmaz.
import os, re, sys, glob

USETS   = [("m10", 3.0, 4.1),      # etiket, U(Mn), U(Ni)   -- taban - 1 eV
           ("p10", 5.0, 6.1)]      #                        -- taban + 1 eV
# 3. set istersen asagidaki satirin # isaretini kaldir (+8 SCF):
# USETS += [("Ni61", 4.0, 6.1)]
STRAINS = [-5]
CFGS    = ["FM", "A", "C", "G"]

def find_base():
    F = {}
    for f in glob.glob("e*_*.in"):
        m = re.match(r"e([+-]?\d+)_(FM|A|C|G)\.in$", os.path.basename(f))
        if m: F[(int(m.group(1)), m.group(2))] = f
    # eps=0 ayri isimlendirilmis: eksik olan her anahtari tek tek doldur
    alt = {"FM": "eps0_cfg0.in", "A": "eps0_cfgA.in", "C": "eps0_cfgC.in", "G": "eps0_cfg1.in"}
    for c, f in alt.items():
        if (0, c) not in F and os.path.exists(f): F[(0, c)] = f
    return F

def fingerprint(fn):
    t = open(fn, errors="ignore").read()
    g = lambda p, d="?": (re.search(p, t).group(1) if re.search(p, t) else d)
    k = re.search(r"K_POINTS[^\n]*\n\s*([\d\s]+)", t)
    return (g(r"nat\s*=\s*(\d+)"), g(r"ntyp\s*=\s*(\d+)"),
            g(r"ecutwfc\s*=\s*([\d.]+)"), g(r"ecutrho\s*=\s*([\d.]+)"),
            " ".join(k.group(1).split()) if k else "?",
            str(len(re.findall(r"^\s*U\s+Mn\d*-3d", t, re.M))),
            str(len(re.findall(r"^\s*U\s+Ni\d*-3d", t, re.M))))

BASE = find_base()
missing = [(e, c) for e in STRAINS for c in CFGS if (e, c) not in BASE]
if missing:
    print("HATA: taban girdiler eksik ->", missing)
    print("bulunanlar:", sorted(BASE.keys())); sys.exit(1)

FP = {k: fingerprint(v) for k, v in BASE.items() if k[0] in STRAINS and k[1] in CFGS}
print("-- taban girdi tutarlilik kontrolu (nat ntyp ecutwfc ecutrho kmesh nUMn nUNi) --")
for k in sorted(FP): print("   e=%+d %-3s %-28s %s" % (k[0], k[1], BASE[k], " ".join(FP[k])))
if len(set(FP.values())) != 1:
    print("*** UYARI: taban girdiler AYNI PARAMETRELERDE DEGIL, %d farkli imza ***" % len(set(FP.values())))
    print("*** Devam etme. U taramasi ancak ayni ayarla anlamli. ***")
    sys.exit(1)
print("   -> 8 girdinin hepsi ayni imzada, devam.")
print()

def patch(txt, uMn, uNi, prefix, outdir):
    out, nMn, nNi = [], 0, 0
    for L in txt.splitlines():
        m = re.match(r"^(\s*U\s+Mn\d*-3d\s+)[\d.eEdD+-]+\s*$", L)
        if m: out.append("%s%.4f" % (m.group(1), uMn)); nMn += 1; continue
        m = re.match(r"^(\s*U\s+Ni\d*-3d\s+)[\d.eEdD+-]+\s*$", L)
        if m: out.append("%s%.4f" % (m.group(1), uNi)); nNi += 1; continue
        if re.match(r"^\s*prefix\s*=", L): out.append("    prefix        = '%s'" % prefix); continue
        if re.match(r"^\s*outdir\s*=", L): out.append("    outdir        = '%s'" % outdir); continue
        out.append(L)
    return "\n".join(out) + "\n", nMn, nNi

jobs = []
print("%-20s %-22s %4s %4s %6s %6s" % ("taban","yeni","nMn","nNi","U_Mn","U_Ni"))
for tag, uMn, uNi in USETS:
    for e in STRAINS:
        for c in CFGS:
            src  = BASE[(e, c)]
            name = "U%s_e%+03d_%s" % (tag, e, c)
            safe = name.replace("+", "p").replace("-", "m")
            txt, nMn, nNi = patch(open(src, errors="ignore").read(), uMn, uNi,
                                  "lnmo_" + safe, "./tmp_" + safe)
            if nMn == 0 or nNi == 0:
                sys.exit("HATA: %s icinde HUBBARD U satiri yok (Mn=%d, Ni=%d)" % (src, nMn, nNi))
            open(name + ".in", "w").write(txt)
            jobs.append(name)
            print("%-20s %-22s %4d %4d %6.2f %6.2f" % (src, name + ".in", nMn, nNi, uMn, uNi))

open("U_jobs.txt", "w").write("\n".join(jobs) + "\n")

SL = ["#!/bin/bash","#SBATCH -J LNMO_U","#SBATCH -p hamsi","#SBATCH -A nparan",
 "#SBATCH -N 1","#SBATCH -n 56","#SBATCH --time=1-00:00:00",
 "#SBATCH --array=1-" + str(len(jobs)) + "%4",
 "#SBATCH --output=slurm_U_%A_%a.out","",
 "export OMP_NUM_THREADS=1","module purge","module load comp/oneapi/2023",
 "module load lib/hdf5/1.14.3-oneapi-2023.0","module load apps/espresso/7.2-oneapi-2023.0","",
 'JOB=$(sed -n "${SLURM_ARRAY_TASK_ID}p" U_jobs.txt)',
 'echo "Basliyor: $JOB $(date)"',
 "mpirun pw.x -npool 4 -in ${JOB}.in > ${JOB}.out","RC=$?",
 'echo "Bitti: $JOB rc=$RC $(date)"',
 "SAFE=$(echo $JOB | tr '+-' 'pm')",
 'if [ $RC -eq 0 ]; then rm -rf ./tmp_${SAFE}; fi',""]
open("run_U.slurm", "w").write("\n".join(SL))

print()
print("%d girdi + U_jobs.txt + run_U.slurm yazildi." % len(jobs))
t = open(jobs[0] + ".in").read()
print("ornek kontrol (%s): nat=%s ntyp=%s calc=%s" % (jobs[0],
      re.search(r"nat\s*=\s*(\d+)", t).group(1), re.search(r"ntyp\s*=\s*(\d+)", t).group(1),
      re.search(r"calculation\s*=\s*'([-\w]+)'", t).group(1)))
nb = len(open(BASE[(STRAINS[0], "FM")], errors="ignore").read().splitlines()); nn = len(t.splitlines())
print("satir sayisi: taban %d -> yeni %d  %s" % (nb, nn,
      "AYNI (sadece deger degisti)" if nb == nn else "*** FARKLI, BAK ***"))
print("maliyet: ~%d dugum-saati, array 1-%d%s4 ile duvar saati ~%.0f saat" % (
      len(jobs), len(jobs), "%", len(jobs) / 4.0))
