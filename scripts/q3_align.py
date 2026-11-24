#!/usr/bin/env python3
# q3_align.py -- Q3 (paper-1 konvansiyonu) <-> J||, Jperp hizalamasi ve fitleri
# Girdi : q3_geo.dat (q3_geo.py cikti) + J_strain.dat (fit_all.py cikti)
# Cikti : q3_align.dat  +  ekrana ~14 satir
import sys, numpy as np

GEO = sys.argv[1] if len(sys.argv) > 1 else "q3_geo.dat"
JDT = sys.argv[2] if len(sys.argv) > 2 else "J_strain.dat"

REF_J  = np.array([11.209,11.524,11.821,12.033,12.152,12.206,12.137,12.046,11.867,11.628,11.291])
REF_JP = np.array([14.140,14.021,13.659,13.090,12.380,11.582,10.722, 9.841, 8.969, 8.117, 7.286])

def numeric_table(fn):
    rows = []
    for L in open(fn, errors="ignore"):
        if L.lstrip().startswith("#"): continue
        p = L.split(); v = []
        for t in p:
            try: v.append(float(t))
            except ValueError: v.append(np.nan)
        if len(v) >= 2 and not np.isnan(v[0]): rows.append(v)
    n = min(len(r) for r in rows)
    return np.array([r[:n] for r in rows], float)

G = numeric_table(GEO)
G = G[np.argsort(G[:, 0])]
eps, Q3, Q3m, Q2 = G[:, 0], G[:, 3], G[:, 7], G[:, 4]
if np.abs(eps).max() <= 1.0: eps = eps * 100.0   # kesir -> yuzde

def pick(T, ref, tol=2e-3):
    T = T[np.argsort(T[:, 0])]
    if T.shape[0] != len(ref): return None
    for c in range(T.shape[1]):
        col = T[:, c]
        if np.all(np.isfinite(col)) and np.abs(col - ref).max() < tol: return col
    return None

src = "J_strain.dat"
try:
    T = numeric_table(JDT)
    Jp, Jq = pick(T, REF_J), pick(T, REF_JP)
    if Jp is None or Jq is None: raise ValueError
except Exception:
    Jp, Jq, src = REF_J.copy(), REF_JP.copy(), "GOMULU REFERANS (durum sec.7g) -- J_strain.dat okunamadi/eslesmedi"

y = np.log(Jp / Jq); Jm = (4*Jp + 2*Jq)/6

def fit(x, z, deg=1):
    c = np.polyfit(x, z, deg); r = np.corrcoef(x, z)[0, 1]
    return c, r, float(np.sqrt(np.mean((z - np.polyval(c, x))**2)))

np.savetxt("q3_align.dat",
    np.column_stack([eps, Q3, Q2, Q3m, Jp, Jq, Jm, y]), fmt="%+9.5f",
    header="eps  Q3_Ni  Q2_Ni  Q3_Mn  J_par  J_perp  J_mean  ln(J_par/J_perp)\n"
           "Q3 = (2/sqrt6)(d_ap-d_eq), paper-1 konvansiyonu; J kaynak: " + src)

print("q3_align.dat yazildi | J kaynagi:", src)
c, r, e = fit(eps, Q3);  print("Q3(eps)            = %+.6f*eps %+.6f      r=%+.6f  Q3=0 @ eps=%+.4f%%" % (c[0], c[1], r, -c[1]/c[0]))
c, r, e = fit(Q3, y);    print("ln(J||/Jp) vs Q3   = %+.4f*Q3 %+.5f       r=%+.6f  RMSE=%.5f" % (c[0], c[1], r, e))
print("                     -> Q3=0'da J||/Jperp = %.4f  (tilt kaynakli artik anizotropi %%%.2f)" % (np.exp(c[1]), 100*(np.exp(c[1])-1)))
c2, r2, e2 = fit(eps, y);print("ln(J||/Jp) vs eps  = %+.6f*eps %+.5f     r=%+.6f  RMSE=%.5f   [istatistiksel esdeger]" % (c2[0], c2[1], r2, e2))
c, r, e = fit(Q3, Jp);   print("J||   vs Q3 lineer : r=%+.4f  -> CIFT kanal (lineer bilesen yok)" % r)
c, r, e = fit(Q3, Jp, 2);print("J||   vs Q3 parabol: %+.3f Q3^2 %+.4f Q3 %+.4f   RMSE=%.4f  tepe Q3=%+.5f" % (c[0], c[1], c[2], e, -c[1]/(2*c[0])))
c, r, e = fit(Q3, Jq);   print("Jperp vs Q3 lineer : %+.4f*Q3 %+.4f   r=%+.6f  RMSE=%.4f  -> TEK kanal" % (c[0], c[1], r, e))
c, r, e = fit(Q3, Jm, 2);print("<J1>  vs Q3 parabol: tepe Q3=%+.5f  (eps=%+.2f%%)" % (-c[1]/(2*c[0]),
      np.polyval(np.polyfit(Q3, eps, 1), -c[1]/(2*c[0]))))
print("Q3 araligi: %+.5f ... %+.5f A   |Q2_Ni|max=%.5f   |Q3_Ni/Q3_Mn| ucta: %.2f / %.2f" % (
      Q3.min(), Q3.max(), np.abs(Q2).max(), abs(Q3[0]/Q3m[0]) if Q3m[0] else np.nan,
      abs(Q3[-1]/Q3m[-1]) if Q3m[-1] else np.nan))
print("paper-1 kontrol: Q3(eps=-5)=%+.5f (P1: +0.170)   Q3(eps=+5)=%+.5f (P1: -0.177)" % (Q3[0], Q3[-1]))
