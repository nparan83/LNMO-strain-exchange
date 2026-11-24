# Fig. 2: robustness - D_par/D_perp(eps) with U band and <J1>(eps) for three U (inputs data/dc_D.dat, data/U_scan.dat)
# Paper: N. Paran, La2NiMnO6 in-plane/out-of-plane exchange under biaxial strain (J. Magn. Magn. Mater., submitted).
# Written as a Google Colab cell; this copy reads the data files of this repository and writes pdf/svg/png to ./figs
# (or $OUT). It prints the values quoted in the manuscript for comparison.
import os
DATA = os.environ.get("DATA", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data"))
# ---- LNMO Fig. 2: dayaniklilik (180 mm) robust ----
BASE = DATA + os.sep

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
    "lines.linewidth":1.0,"legend.frameon":False,"legend.handlelength":1.8,
    "savefig.dpi":600,"savefig.bbox":"tight","savefig.pad_inches":0.02})
CU = ["#92C5DE", "#2166AC", "#053061"]; MK = ["^", "o", "s"]; C_D = "#0072B2"

G = np.loadtxt(BASE + "dc_D.dat")
eps, a_pc, c_pc, DD = G[:,0], G[:,6], G[:,7], G[:,11]
U = np.loadtxt(BASE + "U_scan.dat")
ue, uNi, uJp, uJq, uJm = U[:,0], U[:,2], U[:,3], U[:,4], U[:,6]
ULS = sorted(set(uNi))
LAB = [r"$U_{\rm Ni}=4.1$", r"$U_{\rm Ni}=5.1$ (base)", r"$U_{\rm Ni}=6.1$"]
geo = {e: (a_pc[i]/c_pc[i])**2 for i, e in enumerate(eps)}

fig, (axA, axB) = plt.subplots(1, 2, figsize=(180*MM, 68*MM))
fig.subplots_adjust(wspace=0.30)

# --- (a) D_par/D_perp(eps) + U bandi ---
Dband = {}
for e in (-5., 0., 5.):
    v = [(uJp[(ue==e)&(uNi==u)][0]/uJq[(ue==e)&(uNi==u)][0])*geo[e] for u in ULS]
    Dband[e] = (min(v), max(v))
be = np.array([-5., 0., 5.])
lo = np.interp(eps, be, [Dband[e][0] for e in be])
hi = np.interp(eps, be, [Dband[e][1] for e in be])
axA.fill_between(eps, lo, hi, color=C_D, alpha=0.16, lw=0, zorder=0, label=r"$U\pm1$ eV band")
axA.axhline(1, color="0.75", lw=0.6, zorder=0)
axA.plot(eps, DD, "-", color=C_D, lw=1.1, zorder=2)
axA.plot(eps[eps<0],  DD[eps<0],  "o", ms=4, mfc=C_D,     mec=C_D, zorder=3)
axA.plot(eps[eps>=0], DD[eps>=0], "o", ms=4, mfc="white", mec=C_D, mew=0.9, zorder=3)
axA.plot(-0.559, 1.0, "*", ms=8, mfc="#CC79A7", mec="k", mew=0.4, zorder=4)
axA.annotate(r"$D_\parallel\!=\!D_\perp$" "\n" r"$\varepsilon=-0.56\%$", xy=(-0.559, 1.0),
             xytext=(-4.8, 1.42), fontsize=6.4, linespacing=1.3,
             arrowprops=dict(arrowstyle="-", lw=0.5, color="0.45", shrinkA=1, shrinkB=3))
axA.set_xlabel(r"biaxial strain $\varepsilon$  (%)"); axA.set_ylabel(r"$D_\parallel/D_\perp$")
axA.set_xlim(-5.6, 5.6); axA.set_ylim(0.45, 2.15)
axA.text(0.035, 0.95, "(a)", transform=axA.transAxes, va="top", fontweight="bold")
axA.text(0.97, 0.06, r"$2.80^{+0.50}_{-0.31}\times$ total", transform=axA.transAxes,
         ha="right", fontsize=6.8)
axA.legend(loc="upper left", bbox_to_anchor=(0.03, 0.88), borderaxespad=0)

# --- (b) <J1>(eps) uc U icin; tepe konumu sabit ---
x = np.linspace(-5.6, 5.6, 300); PK = []
for k, u in enumerate(ULS):
    m = uNi == u; ee, yy = ue[m], uJm[m]
    o = np.argsort(ee); ee, yy = ee[o], yy[o]
    c = np.polyfit(ee, yy, 2); pk = -c[1]/(2*c[0]); PK.append(pk)
    axB.plot(x, np.polyval(c, x), "-", color=CU[k], lw=1.0, zorder=1)
    axB.plot(ee, yy, MK[k], ms=4.2, mfc=CU[k], mec="k", mew=0.4, ls="none",
             zorder=3, label=LAB[k])
    axB.plot([pk, pk], [7.9, np.polyval(c, pk)], ":", color=CU[k], lw=0.8, zorder=2)
    axB.plot(pk, np.polyval(c, pk), "v", ms=3.4, mfc=CU[k], mec="k", mew=0.35, zorder=4)
axB.set_xlabel(r"biaxial strain $\varepsilon$  (%)"); axB.set_ylabel(r"$\langle J_1\rangle$  (meV)")
axB.set_xlim(-5.6, 5.6); axB.set_ylim(7.9, 14.9)
axB.text(0.035, 0.95, "(b)", transform=axB.transAxes, va="top", fontweight="bold")
axB.legend(loc="lower center", borderaxespad=0.6, labelspacing=0.35)
axB.annotate(r"$\varepsilon_{\rm peak} = %+.2f\ /\ %+.2f\ /\ %+.2f\ \%%$" % tuple(PK) + "\n"
             r"position invariant ($\pm0.03\%$)" "\n" r"magnitude varies $\pm12\%$",
             xy=(-3.0, 13.95), xytext=(-0.7, 14.73), fontsize=5.6, linespacing=1.3,
             va="top", ha="left",
             arrowprops=dict(arrowstyle="->", lw=0.5, color="0.45", shrinkA=3, shrinkB=2))

for _ext, _kw in (("pdf", {}), ("svg", {}), ("png", {"dpi": 600})):
    _p = os.path.join(OUT, "fig_robust." + _ext); fig.savefig(_p, **_kw); print("yazildi:", _p)
plt.show()
print("tepe eps:", " / ".join("%+.3f%%" % p for p in PK))
print("D band:", {k: (round(v[0],3), round(v[1],3)) for k, v in Dband.items()})

from PIL import Image
_p = os.path.join(OUT, "fig_robust.png")
w, h = Image.open(_p).size
print(os.path.basename(_p), w, "x", h, "px  ->  %.0f dpi @ 180 mm" % (w/(180/25.4)))
