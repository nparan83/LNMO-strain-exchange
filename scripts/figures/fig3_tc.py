# Fig. 3: Binder cumulant crossing and T_C(eps) (inputs data/mc/mc_eps*.dat)
# Paper: N. Paran, La2NiMnO6 in-plane/out-of-plane exchange under biaxial strain (J. Magn. Magn. Mater., submitted).
# Written as a Google Colab cell; this copy reads the data files of this repository and writes pdf/svg/png to ./figs
# (or $OUT). It prints the values quoted in the manuscript for comparison.
import os
DATA = os.environ.get("DATA", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data"))
# ---- LNMO Fig. 3: T_C ve Binder (180 mm) ----
BASE = DATA + os.sep

import glob, re
import numpy as np, matplotlib as mpl, matplotlib.pyplot as plt

import matplotlib as mpl
mpl.rcParams['figure.dpi']  = 200                # SADECE ekran; dosyayi etkilemez
mpl.rcParams['savefig.dpi'] = 600

import os
OUT = os.environ.get("OUT", "figs")
os.makedirs(OUT, exist_ok=True)

MM = 1/25.4
mpl.rcParams.update({
    "pdf.fonttype":42, "ps.fonttype":42, "svg.fonttype":"path",
    "font.family":"sans-serif","font.sans-serif":["DejaVu Sans","Arial"],
    "font.size":8,"axes.labelsize":8,"legend.fontsize":6.6,
    "xtick.labelsize":7,"ytick.labelsize":7,"axes.linewidth":0.7,
    "xtick.direction":"in","ytick.direction":"in","xtick.top":True,"ytick.right":True,
    "xtick.major.width":0.7,"ytick.major.width":0.7,"xtick.minor.visible":True,
    "ytick.minor.visible":True,"xtick.minor.width":0.5,"ytick.minor.width":0.5,
    "lines.linewidth":1.0,"legend.frameon":False,"legend.handlelength":1.6,
    "savefig.dpi":600,"savefig.bbox":"tight","savefig.pad_inches":0.02})
CL = ["#BDD7E7", "#6BAED6", "#2171B5", "#08306B"]   # L artan: acik -> koyu
ML = ["o", "s", "^", "D"]
C_TC, C_EXP = "#0072B2", "#009E73"
PANEL_EPS = 0            # (a) panelinde gosterilecek strain

def load(fn):
    D = np.loadtxt(fn)
    e = int(re.search(r"eps=([+-]?\d+)", open(fn).readline()).group(1))   # ADDAN DEGIL
    return e, D

def crossings(D):
    L, T, U4 = D[:,0], D[:,1], D[:,3]
    Ls = sorted(set(L.astype(int)))
    cur = {l: (np.sort(T[L==l]), U4[L==l][np.argsort(T[L==l])]) for l in Ls}
    cr = []
    for i in range(len(Ls)):
        for j in range(i+1, len(Ls)):
            t1,u1 = cur[Ls[i]]; t2,u2 = cur[Ls[j]]
            tg = np.linspace(max(t1.min(),t2.min()), min(t1.max(),t2.max()), 4000)
            d = np.interp(tg,t1,u1) - np.interp(tg,t2,u2)
            s = np.where(np.diff(np.sign(d)))[0]
            if len(s):
                k = s[0]; tc = tg[k] - d[k]*(tg[k+1]-tg[k])/(d[k+1]-d[k])
                cr.append((tc, np.interp(tc,t1,u1)))
    a = np.array(cr)
    return a[:,0].mean(), a[:,0].std(ddof=1), a[:,1].mean(), cur, Ls

REC = {}
for fn in sorted(glob.glob(os.path.join(DATA, "mc", "mc_eps*.dat"))):
    e, D = load(fn); REC[e] = crossings(D)
E = np.array(sorted(REC)); TC = np.array([REC[e][0] for e in E])
SD = np.array([REC[e][1] for e in E]); U4S = np.array([REC[e][2] for e in E])

REF = {-5:304.3,-4:309.9,-3:312.3,-2:310.2,-1:306.4,0:299.9,1:294.0,2:283.5,3:271.6,4:259.4,5:246.6}
bad = [(int(e), round(t,2), REF[int(e)]) for e, t in zip(E, TC) if abs(t-REF[int(e)]) > 2.0]
print("sec.7h karsilastirmasi (tol 2 K):", "11/11 DOGRULANDI" if not bad else "*** SAPMA %s" % bad)
print("U4* = %.4f +- %.4f   (11 noktada sabit; 2/3 ve 4/9 limitleri arasinda)" % (U4S.mean(), U4S.std(ddof=1)))

fig, (axA, axB) = plt.subplots(1, 2, figsize=(180*MM, 70*MM))
fig.subplots_adjust(wspace=0.30)

# ---- (a) Binder kesisimi ----
tcm, tcs, u4m, cur, Ls = REC[PANEL_EPS]
for k, l in enumerate(Ls):
    t, u = cur[l]
    axA.plot(t, u, "-", color=CL[k], lw=0.9, zorder=1)
    axA.plot(t, u, ML[k], ms=3.0, mfc=CL[k], mec=CL[k], ls="none", zorder=2, label=r"$L=%d$" % l)
axA.axhline(2/3, color="0.75", lw=0.6, ls="--", zorder=0)
axA.axhline(4/9, color="0.75", lw=0.6, ls="--", zorder=0)
axA.text(t.max(), 2/3+0.004, r"$2/3$  ($T\!\to\!0$)", fontsize=6.0, color="0.45",
         ha="right", va="bottom")
axA.text(t.max(), 4/9-0.004, r"$4/9$  ($T\!\to\!\infty$)", fontsize=6.0, color="0.45",
         ha="right", va="top")
axA.axvline(tcm, color=C_TC, lw=0.7, ls=":", zorder=0)
axA.set_xlabel(r"$T$  (K)"); axA.set_ylabel(r"Binder cumulant  $U_4$")
axA.set_ylim(0.418, 0.692)
axA.text(0.02, 0.98, "(a)", transform=axA.transAxes, va="top", fontweight="bold")
axA.text(0.035, 0.868, r"$\varepsilon=%+d\%%$" % PANEL_EPS, transform=axA.transAxes,
         va="top", fontsize=7)
axA.legend(loc="upper right", bbox_to_anchor=(0.985, 0.90), labelspacing=0.28,
           ncol=2, columnspacing=0.9, handletextpad=0.5)

axI = axA.inset_axes([0.115, 0.115, 0.42, 0.315])
for k, l in enumerate(Ls):
    t, u = cur[l]
    axI.plot(t, u, "-", color=CL[k], lw=0.9)
    axI.plot(t, u, ML[k], ms=2.6, mfc=CL[k], mec=CL[k], ls="none")
axI.plot(tcm, u4m, "*", ms=8, mfc="#CC79A7", mec="k", mew=0.4, zorder=5)
axI.set_xlim(tcm-6, tcm+6); axI.set_ylim(u4m-0.013, u4m+0.013)
axI.tick_params(labelsize=5.4, length=2.2, pad=1.5)
for sp in axI.spines.values(): sp.set_linewidth(0.6)
axI.text(0.045, 0.055, r"$T_{\rm C}=%.1f\pm%.1f$ K" % (tcm, tcs) + "\n"
         + r"$U_4^{*}=%.4f$" % u4m, transform=axI.transAxes, fontsize=5.5,
         va="bottom", ha="left", linespacing=1.35)

# ---- (b) T_C(eps) ----
axB.axhspan(273, 280, color=C_EXP, alpha=0.18, lw=0, zorder=0)
axB.text(5.28, 276.5, "exp.", fontsize=6.2, color="#00654A", ha="right", va="center")
# MODEL SECIMI: parabol artiklari yapisal (RMSE 1.14 K, MC hatasi 0.65 K);
# kubik RMSE 0.34 K, artiklar gurultu seviyesinde -> kubik.
x = np.linspace(-5.5, 5.5, 300); c = np.polyfit(E, TC, 3)
for dg in (2, 3, 4):
    cc = np.polyfit(E, TC, dg); rr = TC - np.polyval(cc, E)
    print("  derece %d: RMSE %.2f K  max|artik| %.2f K" % (dg, np.sqrt((rr**2).mean()), np.abs(rr).max()))
axB.plot(x, np.polyval(c, x), "-", color=C_TC, lw=1.0, zorder=1)
axB.errorbar(E[E<0], TC[E<0], yerr=SD[E<0], fmt="o", ms=4, mfc=C_TC, mec=C_TC,
             ecolor=C_TC, elinewidth=0.8, capsize=1.8, capthick=0.8, zorder=3)
axB.errorbar(E[E>=0], TC[E>=0], yerr=SD[E>=0], fmt="o", ms=4, mfc="white", mec=C_TC,
             mew=0.9, ecolor=C_TC, elinewidth=0.8, capsize=1.8, capthick=0.8, zorder=3)
dc = np.polyder(c); pk = [r.real for r in np.roots(dc) if abs(r.imag) < 1e-9 and -5 < r.real < 5]
pk = pk[int(np.argmax(np.polyval(c, pk)))]
axB.plot([pk, pk], [238, np.polyval(c, pk)], ":", color="#CC79A7", lw=0.9, zorder=2)
axB.plot(pk, np.polyval(c, pk), "v", ms=4.2, mfc="#CC79A7", mec="k", mew=0.4, zorder=4)
axB.annotate(r"peak  $\varepsilon=%+.2f\%%$" % pk, xy=(pk-0.12, np.polyval(c, pk)),
             xytext=(-0.30, 315), fontsize=6.4, ha="left",
             arrowprops=dict(arrowstyle="-", lw=0.5, color="0.45", shrinkA=1, shrinkB=3))
axB.set_xlabel(r"biaxial strain $\varepsilon$  (%)"); axB.set_ylabel(r"$T_{\rm C}$  (K)")
axB.set_xlim(-5.5, 5.5); axB.set_ylim(238, 321)
axB.text(0.02, 0.98, "(b)", transform=axB.transAxes, va="top", fontweight="bold")
axB.text(0.035, 0.10, r"$\Delta T_{\rm C}=66$ K  (22%)" "\n"
         r"$\pm12\%$ systematic ($U$; Fig. 2b)", transform=axB.transAxes,
         fontsize=6.2, linespacing=1.35)

for _ext, _kw in (("pdf", {}), ("svg", {}), ("png", {"dpi": 600})):
    _p = os.path.join(OUT, "fig_tc." + _ext); fig.savefig(_p, **_kw); print("yazildi:", _p)
plt.show()
print("T_C tepe (kubik fit): %+.3f%%   T_C araligi %.1f - %.1f K" % (pk, TC.min(), TC.max()))

from PIL import Image
_p = os.path.join(OUT, "fig_tc.png")
w, h = Image.open(_p).size
print(os.path.basename(_p), w, "x", h, "px  ->  %.0f dpi @ 180 mm" % (w/(180/25.4)))
