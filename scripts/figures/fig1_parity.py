# Fig. 1: parity of J_par and J_perp in Q3 (input data/q3_align.dat)
# Paper: N. Paran, La2NiMnO6 in-plane/out-of-plane exchange under biaxial strain (J. Magn. Magn. Mater., submitted).
# Written as a Google Colab cell; this copy reads the data files of this repository and writes pdf/svg/png to ./figs
# (or $OUT). It prints the values quoted in the manuscript for comparison.
import os
DATA = os.environ.get("DATA", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data"))
# ---- LNMO Fig. 1: parite ayrismasi (180 mm) ----
PATH = os.path.join(DATA, "q3_align.dat")

import numpy as np, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

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
    "font.size":8,"axes.labelsize":8,"axes.titlesize":8,"legend.fontsize":7,
    "xtick.labelsize":7,"ytick.labelsize":7,"axes.linewidth":0.7,
    "xtick.direction":"in","ytick.direction":"in","xtick.top":True,"ytick.right":True,
    "xtick.major.width":0.7,"ytick.major.width":0.7,"xtick.minor.visible":True,
    "ytick.minor.visible":True,"xtick.minor.width":0.5,"ytick.minor.width":0.5,
    "lines.linewidth":1.0,"legend.frameon":False,"legend.handlelength":1.8,
    "savefig.dpi":600,"savefig.bbox":"tight","savefig.pad_inches":0.02,
})
C_PAR, C_PERP, C_FIT = "#0072B2", "#D55E00", "#333333"   # Okabe-Ito

D = np.loadtxt(PATH)
eps, Q3, Jp, Jq, ln = D[:,0], D[:,1], D[:,4], D[:,5], D[:,7]
if np.abs(eps).max() <= 1.0: eps = eps*100.0
comp = eps < 0
q  = np.linspace(Q3.min()-0.008, Q3.max()+0.008, 300)
cP = np.polyfit(Q3, Jp, 2); cQ = np.polyfit(Q3, Jq, 1)
c1 = np.polyfit(Q3, ln, 1); c2 = np.polyfit(Q3, ln, 2)
r1 = ln - np.polyval(c1, Q3); r2 = ln - np.polyval(c2, Q3)
ke = np.polyfit(Q3, eps, 1)

fig = plt.figure(figsize=(180*MM, 72*MM))
gs  = GridSpec(2, 2, height_ratios=[3.0, 1.0], hspace=0.10, wspace=0.30, figure=fig)
axA = fig.add_subplot(gs[:, 0]); axB = fig.add_subplot(gs[0, 1])
axR = fig.add_subplot(gs[1, 1], sharex=axB)

def mk(ax, x, y, c, m="o"):
    ax.plot(x[comp],  y[comp],  m, ms=4, mfc=c,       mec=c, mew=0.8, ls="none")
    ax.plot(x[~comp], y[~comp], m, ms=4, mfc="white", mec=c, mew=0.9, ls="none")

axA.axvline(0, color="0.75", lw=0.6, zorder=0)
axA.plot(q, np.polyval(cP, q), "-",  color=C_PAR,  lw=1.0, zorder=1)
axA.plot(q, np.polyval(cQ, q), "--", color=C_PERP, lw=1.0, zorder=1)
mk(axA, Q3, Jp, C_PAR, "o"); mk(axA, Q3, Jq, C_PERP, "s")
axA.set_xlabel(r"$Q_3$ (NiO$_6$)  ($\mathrm{\AA}$)"); axA.set_ylabel(r"$J_1$  (meV)")
axA.set_xlim(q[0], q[-1]); axA.set_ylim(6.3, 15.0)
h = [plt.Line2D([],[], color=C_PAR,  marker="o", ms=4, lw=1.0,
                label=r"$J_\parallel$  (even: $\propto Q_3^{2}$)"),
     plt.Line2D([],[], color=C_PERP, marker="s", ms=4, lw=1.0, ls="--",
                label=r"$J_\perp$  (odd: $\propto Q_3$)"),
     plt.Line2D([],[], color="0.35", marker="o", ms=4, ls="none", mfc="0.35",
                label=r"compressive  ($\varepsilon<0$)"),
     plt.Line2D([],[], color="0.35", marker="o", ms=4, ls="none", mfc="white", mew=0.9,
                label=r"tensile  ($\varepsilon>0$)")]
axA.legend(handles=h, loc="lower right", borderaxespad=0.5, labelspacing=0.35)
axA.text(0.03, 0.95, "(a)", transform=axA.transAxes, va="top", fontweight="bold")
axA.text(0.004, 14.55, r"$Q_3=0$", fontsize=6.5, color="0.45")
f2e = (lambda z: np.polyval(ke, z), lambda e: (e - ke[1]) / ke[0])
axA.secondary_xaxis("top", functions=f2e).set_xlabel(r"biaxial strain $\varepsilon$  (%)", labelpad=4)

axB.axhline(0, color="0.75", lw=0.6, zorder=0); axB.axvline(0, color="0.75", lw=0.6, zorder=0)
axB.plot(q, np.polyval(c2, q), "-", color=C_FIT, lw=1.0, zorder=1)
axB.plot(q, np.polyval(c1, q), ":", color="0.55", lw=0.9, zorder=1)
mk(axB, Q3, ln, "#009E73", "D")
axB.set_ylabel(r"$\ln(J_\parallel/J_\perp)$", labelpad=1)
axB.set_ylim(-0.31, 0.50); plt.setp(axB.get_xticklabels(), visible=False)
axB.text(0.1, 0.93, "(b)", transform=axB.transAxes, va="top", fontweight="bold")
axB.plot([], [], "-", color=C_FIT, lw=1.0, label="quadratic")
axB.plot([], [], ":", color="0.55", lw=0.9, label="linear")
axB.legend(loc="lower left", borderaxespad=0.5, labelspacing=0.3)
axB.plot(0, c2[2], "*", ms=7, mfc="#CC79A7", mec="k", mew=0.4, zorder=5)
axB.annotate(r"$J_\parallel/J_\perp = %.3f$" % np.exp(c2[2]) + "\n" + r"at $Q_3 = 0$",
             xy=(0.006, c2[2]+0.012), xytext=(0.045, 0.30), fontsize=6.3, va="center",
             linespacing=1.3, arrowprops=dict(arrowstyle="-", lw=0.5, color="0.45",
                                              shrinkA=1, shrinkB=2))
axB.secondary_xaxis("top", functions=f2e).set_xlabel(r"$\varepsilon$  (%)", labelpad=4)

axR.axhline(0, color="0.75", lw=0.6)
axR.plot(Q3, r1, ":", color="0.55", lw=0.9); axR.plot(Q3, r1, "o", ms=2.6, mfc="0.55", mec="0.55")
axR.plot(Q3, r2, "-", color=C_FIT, lw=0.9);  axR.plot(Q3, r2, "o", ms=2.6, mfc=C_FIT, mec=C_FIT)
axR.set_xlabel(r"$Q_3$ (NiO$_6$)  ($\mathrm{\AA}$)"); axR.set_ylabel("residual", labelpad=1)
axR.set_ylim(-0.055, 0.055); axR.set_yticks([-0.04, 0, 0.04]); axR.set_xlim(q[0], q[-1])

for _ext, _kw in (("pdf", {}), ("svg", {}), ("png", {"dpi": 600})):
    _p = os.path.join(OUT, "fig_parity." + _ext); fig.savefig(_p, **_kw); print("yazildi:", _p)
plt.show()
print("J_par  = %+.4f Q3^2 %+.4f Q3 %+.4f" % tuple(cP))
print("J_perp = %+.4f Q3 %+.4f   r = %.5f" % (*cQ, np.corrcoef(Q3, Jq)[0,1]))
print("ln lin : %+.4f Q3 %+.5f   RMSE %.5f" % (*c1, np.sqrt((r1**2).mean())))
print("ln quad: %+.4f Q3^2 %+.4f Q3 %+.5f   RMSE %.5f" % (*c2, np.sqrt((r2**2).mean())))
print("Q3=0 -> J_par/J_perp = %.4f" % np.exp(c2[2]))

from PIL import Image
_p = os.path.join(OUT, "fig_parity.png")
w, h = Image.open(_p).size
print(os.path.basename(_p), w, "x", h, "px  ->  %.0f dpi @ 180 mm" % (w/(180/25.4)))
