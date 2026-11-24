# Fig. 4: static versus relaxed energy mapping and E_rel (input data/Erel.dat)
# Paper: N. Paran, La2NiMnO6 in-plane/out-of-plane exchange under biaxial strain (J. Magn. Magn. Mater., submitted).
# Written as a Google Colab cell; this copy reads the data files of this repository and writes pdf/svg/png to ./figs
# (or $OUT). It prints the values quoted in the manuscript for comparison.
import os
DATA = os.environ.get("DATA", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data"))
# ---- LNMO Fig. 4: statik vs gevsetilmis + E_rel (180 mm) ----
PATH = os.path.join(DATA, "Erel.dat")

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
    "font.size":8,"axes.labelsize":8,"legend.fontsize":6.4,
    "xtick.labelsize":7,"ytick.labelsize":7,"axes.linewidth":0.7,
    "xtick.direction":"in","ytick.direction":"in","xtick.top":True,"ytick.right":True,
    "xtick.major.width":0.7,"ytick.major.width":0.7,"xtick.minor.visible":True,
    "ytick.minor.visible":True,"xtick.minor.width":0.5,"ytick.minor.width":0.5,
    "lines.linewidth":1.0,"legend.frameon":False,"legend.handlelength":1.8,
    "savefig.dpi":600,"savefig.bbox":"tight","savefig.pad_inches":0.02})
C_C, C_G, C_V = "#0072B2", "#D55E00", "#009E73"

D = np.loadtxt(PATH)
eps, V_FM = D[:,0], D[:,2]
dE_AF, Cst, relC, dV_AF = D[:,5], D[:,6], D[:,7], D[:,8]
dE_Fi, Gst, relG, dV_Fi = D[:,11], D[:,12], D[:,13], D[:,14]
g = np.isfinite(dE_Fi)

fig, (axA, axB) = plt.subplots(1, 2, figsize=(180*MM, 70*MM))
fig.subplots_adjust(wspace=0.34)

# ---- (a) statik vs gevsetilmis ----
axA.fill_between(eps, dE_AF, Cst, color=C_C, alpha=0.16, lw=0, zorder=0)
axA.fill_between(eps[g], dE_Fi[g], Gst[g], color=C_G, alpha=0.16, lw=0, zorder=0)
axA.plot(eps, Cst, "--", color=C_C, lw=1.0, zorder=1)
axA.plot(eps, dE_AF, "-", color=C_C, lw=1.0, zorder=2)
axA.plot(eps, dE_AF, "o", ms=3.6, mfc=C_C, mec=C_C, ls="none", zorder=3)
axA.plot(eps, Gst, "--", color=C_G, lw=1.0, zorder=1)
axA.plot(eps[g], dE_Fi[g], "-", color=C_G, lw=1.0, zorder=2)
axA.plot(eps[g], dE_Fi[g], "s", ms=3.8, mfc=C_G, mec=C_G, ls="none", zorder=3)
axA.set_xlabel(r"biaxial strain $\varepsilon$  (%)")
axA.set_ylabel(r"$\Delta E$ vs FM  (meV/f.u.)")
axA.set_xlim(-5.5, 5.5); axA.set_ylim(100, 243)
axA.text(0.025, 0.98, "(a)", transform=axA.transAxes, va="top", fontweight="bold")
axA.text(1, 227, r"G  (FiM,  $M_{\rm tot}=2\,\mu_B$)", fontsize=6.6, color=C_G, ha="center")
axA.text(0, 120, r"C  (AFM,  $M_{\rm tot}=0$)", fontsize=6.6, color=C_C, ha="center")
h = [plt.Line2D([],[], color="0.35", ls="--", lw=1.0, label="static (fixed FM geometry)"),
     plt.Line2D([],[], color="0.35", ls="-",  lw=1.0, marker="o", ms=3.6,
                label="relaxed (vc-relax)"),
     plt.Rectangle((0,0),1,1, fc="0.55", alpha=0.25, ec="none",
                   label=r"$E_{\rm rel}$  (magnetostructural)")]
axA.legend(handles=h, loc="upper right", bbox_to_anchor=(0.7, 0.550),
           borderaxespad=0.3, labelspacing=0.32)

# ---- (b) E_rel ve Dc/c ----
axB.plot(eps, relC, "-", color=C_C, lw=1.0, zorder=2)
axB.plot(eps, relC, "o", ms=4, mfc=C_C, mec=C_C, ls="none", zorder=3)
axB.plot(eps[g], relG[g], "--", color=C_G, lw=1.0, zorder=2)
axB.plot(eps[g], relG[g], "s", ms=4.2, mfc=C_G, mec=C_G, ls="none", zorder=3)
axB.set_xlabel(r"biaxial strain $\varepsilon$  (%)")
axB.set_ylabel(r"$E_{\rm rel}=\Delta E_{\rm static}-\Delta E_{\rm relaxed}$  (meV/f.u.)")
axB.set_xlim(-5.5, 5.5); axB.set_ylim(9.0, 28.5)
axB.text(0.025, 0.98, "(b)", transform=axB.transAxes, va="top", fontweight="bold")
hb = [plt.Line2D([],[], color=C_C, marker="o", ms=4, lw=1.0, label=r"$E_{\rm rel}$(C),  11 pts"),
      plt.Line2D([],[], color=C_G, marker="s", ms=4.2, lw=1.0, ls="--", label=r"$E_{\rm rel}$(G),  3 pts"),
      plt.Line2D([],[], color=C_V, marker="^", ms=3.4, lw=0.9, ls=":", mfc="white", mew=0.9,
                 label=r"$\Delta c_\perp/c_\perp$ (C)"),
      plt.Line2D([],[], color=C_V, marker="v", ms=3.6, lw=0, ls="none",
                 label=r"$\Delta c_\perp/c_\perp$ (G)")]
axB.legend(handles=hb, loc="center left", bbox_to_anchor=(0.57, 0.86),
           labelspacing=0.32, handletextpad=0.6)

axT = axB.twinx()
axT.plot(eps, 100*dV_AF/V_FM, "^", ms=3.4, mfc="white", mec=C_V, mew=0.9, ls="none", zorder=3)
axT.plot(eps, 100*dV_AF/V_FM, ":", color=C_V, lw=0.9, zorder=2)
axT.plot(eps[g], 100*dV_Fi[g]/V_FM[g], "v", ms=3.6, mfc=C_V, mec=C_V, ls="none", zorder=3)
axT.set_ylabel(r"$\Delta c_\perp/c_\perp$  (%)", color=C_V)
axT.tick_params(axis="y", colors=C_V, direction="in"); axT.spines["right"].set_color(C_V)
axT.set_ylim(-0.03, 0.62); axT.minorticks_on()
axB.tick_params(right=False)

for _ext, _kw in (("pdf", {}), ("svg", {}), ("png", {"dpi": 600})):
    _p = os.path.join(OUT, "fig_erel." + _ext); fig.savefig(_p, **_kw); print("yazildi:", _p)
plt.show()
print("E_rel(C) = %.2f +- %.2f meV/f.u.  |  E_rel(G) = %s" % (
      relC.mean(), relC.std(ddof=1), " / ".join("%.2f" % v for v in relG[g])))
print("kirilan bag basina: C %.2f  G %.2f meV" % (relC.mean()/8, relG[g].mean()/12))

from PIL import Image
_p = os.path.join(OUT, "fig_erel.png")
w, h = Image.open(_p).size
print(os.path.basename(_p), w, "x", h, "px  ->  %.0f dpi @ 180 mm" % (w/(180/25.4)))
