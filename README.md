# LNMO-strain-exchange

Data and code for

> N. Paran, *Opposite parity of the in-plane and out-of-plane Ni–O–Mn exchange displaces the Curie maximum of biaxially strained La₂NiMnO₆ toward compression*, submitted to Journal of Magnetism and Magnetic Materials (2026).

Companion data set (relaxed structures, band structure and the C-type/G-type vc-relax runs used here): [LNMO-strain-parity](https://github.com/nparan83/LNMO-strain-parity), doi:[10.5281/zenodo.24353930](https://doi.org/10.5281/zenodo.24353930).

## Calculations

Quantum ESPRESSO 7.2, PBE+U (Dudarev, ortho-atomic projectors), U(Mn) = 4.0 eV, U(Ni) = 5.1 eV, PSlibrary PAW datasets La.pbe-spfn-kjpaw_psl.1.0.0, Mn.pbe-spn-kjpaw_psl.0.3.1, Ni.pbe-n-kjpaw_psl.0.1, O.pbe-n-kjpaw_psl.1.0.0, cut-offs 90/900 Ry, Marzari–Vanderbilt smearing 0.01 Ry, conv_thr 10⁻⁸ Ry, 40-atom √2×√2×1 supercell of the relaxed P2₁/n cell, 4×4×4 k-mesh. All spin configurations are static SCF runs at the relaxed ferromagnetic geometry of each strain state. The pseudopotential files are not redistributed; they are available from PSlibrary.

Exchange model: H = −Σ⟨ij⟩ J_ij **S**_i·**S**_j, each bond counted once, |S(Ni)| = 1, |S(Mn)| = 3/2, J > 0 ferromagnetic. The 24 first-neighbour Ni–O–Mn bonds of the supercell split into 16 in-plane (J∥) and 8 out-of-plane (J⊥) bonds. The four configurations form a 2² factorial design; with E in meV per 40-atom cell,

| configuration | sign rule | E |
|---|---|---|
| FM | +1 | E₀ − 24 J∥ − 12 J⊥ |
| A | (−1)^m₃ | E₀ − 24 J∥ + 12 J⊥ |
| C | (−1)^(m₁+m₂) | E₀ + 24 J∥ − 12 J⊥ |
| G | (−1)^(m₁+m₂+m₃) | E₀ + 24 J∥ + 12 J⊥ |

Monte Carlo: classical anisotropic Heisenberg model on the simple-cubic Ni/Mn rock-salt lattice (x, y bonds J∥, z bonds J⊥), Metropolis, L = 8, 12, 16, 24, 2×10⁴ equilibration and 1.2×10⁵ measurement sweeps, Binder-cumulant crossings.

## Layout

| Folder | Content |
|---|---|
| `inputs/strain_2x2/`, `outputs/strain_2x2/` | 44 SCF runs: 11 strain states × {FM, A, C, G} |
| `inputs/U_scan/`, `outputs/U_scan/` | 24 SCF runs: U(Mn)/U(Ni) = 3.0/4.1 (`Um10`) and 5.0/6.1 eV (`Up10`) at ε = −5, 0, +5 % |
| `inputs/pilot_isotropic/`, `outputs/pilot_isotropic/` | isotropic pilot at ε = 0 (configurations cfg2–cfg5; cfg0 and cfg1 are in `strain_2x2`) |
| `inputs/kpoint_test/`, `outputs/kpoint_test/` | FM and G at ε = 0 on a 6×6×6 mesh |
| `inputs/static_G_eps0/`, `outputs/static_G_eps0/` | 20-atom G-type SCF at the FM geometry (cross-check, 215.97 meV/f.u.) |
| `data/` | derived data files (columns are given in the header line of each file) |
| `data/mc/` | Monte Carlo output `mc_eps<ε>.dat`: L, T, ⟨m⟩, U₄, error (ε in the header line) |
| `scripts/` | generation, fitting, Monte Carlo and table scripts (Python) |
| `scripts/figures/` | Fig. 1–4 |
| `slurm/` | SLURM job scripts used on the TRUBA cluster |
| `MANIFEST.txt` | size and path of every file of the original package |

Output files are gzip-compressed.

**File names.** `e-05_FM.in` … `e+05_G.in` are ε = −5 … +5 % (in percent; ε > 0 is tensile). At ε = 0 the runs are named `eps0_cfg0` (FM), `eps0_cfgA` (A), `eps0_cfgC` (C) and `eps0_cfg1` (G). In the 40-atom inputs the species are grouped: types 3–6 are Mn1–Mn4 and types 7–10 are Ni1–Ni4.

## From files to results

| Manuscript item | Data file | Produced by |
|---|---|---|
| Sec. 3.1–3.2, Table 2: J∥, J⊥, ⟨J₁⟩, residual, moments | `data/J_strain.dat` | `scripts/fit_all.py` (reads `outputs/strain_2x2`) |
| Sec. 2.2, Table 1: Ni–O and Mn–O distances, θ∥, θ⊥, Q₃, Q₂ | `data/q3_geo.dat`, `data/geo_theta.dat` | `scripts/q3_geo.py` (uses `q3lib.py`; reads the relaxed FM structures of the companion data set), `scripts/geo2.py` |
| Sec. 3.2, Fig. 1: J and ln(J∥/J⊥) against Q₃ | `data/q3_align.dat` | `scripts/q3_align.py`; figure `scripts/figures/fig1_parity.py` |
| Sec. 2.2, Table 1 and Sec. 3.3, Fig. 2(a), Table 2: lattice parameters, β, a_pc, c⊥, D∥/D⊥ | `data/dc_D.dat` | `scripts/getc.py`; figure `scripts/figures/fig2_robust.py` |
| Sec. 3.4, Table 3, Fig. 3: T_C(ε), T_MFA | `data/Tc_strain.dat`, `data/mc/mc_eps*.dat` | `scripts/mc2.py` (`slurm/run_mc.slurm`), `scripts/tc.py`; figure `scripts/figures/fig3_tc.py` |
| Sec. 3.6, Fig. 4, Table 1 (volume): E_rel(ε), FM elastic curve, ΔV, Δc⊥ | `data/Erel.dat` | `scripts/mkerel.py` (reads the C-type and G-type vc-relax outputs of the companion data set); figure `scripts/figures/fig4_erel.py` |
| Sec. 3.7, Table 4, Fig. 2(b): Hubbard-U scan | `data/U_scan.dat` | `scripts/mkU.py` (inputs), `scripts/fit_U.py` |
| Tables 1–4 (Appendices B and C) | `data/tables.md`, `data/tables.tex` | `scripts/mktab.py` (reads the six files above; prints 12 reference checks) |
| Sec. 2.3: supercell and configurations | `inputs/strain_2x2/` | `scripts/mkstrain.py`, `scripts/mkAC2.py`; check `scripts/afm_id2.py` |
| Sec. 2.3: isotropic pilot, k-mesh test | `outputs/pilot_isotropic/`, `outputs/kpoint_test/`, `data/design.npz` | `scripts/fit_J.py`, `scripts/fit_JJ.py`, `scripts/fit_pm5.py` |
| Table 5: computational cost | wall times in the `.out` files | — |

The substrate mismatches of Sec. 4 are taken from the companion article (Table 5 there, a_pc = 3.935 Å).

## Reproducing the derived data

```
bash scripts/reproduce_all.sh                       # from the outputs in this repository
bash scripts/reproduce_all.sh ../LNMO-strain-parity # also rebuilds q3_geo.dat, dc_D.dat and Erel.dat
```

The script unpacks the outputs into `./work_reproduce`, runs `fit_all.py`, `q3_align.py`, `tc.py`, `fit_U.py` and `mktab.py` (and, with the companion data set, `q3_geo.py`, `getc.py`, `mkerel.py`), and compares every regenerated file with `data/`. Both variants reproduce all nine files exactly (checked with Python 3.13, NumPy 2). The Monte Carlo itself (`mc2.py`) is not rerun; its outputs are `data/mc/`. The figure scripts write pdf/svg/png to `./figs` and print the values quoted in the text.

Two notes on the Monte Carlo analysis. `tc.py` takes T_C from the crossing of the L = 16 and L = 24 cumulants (Table 3). `fig3_tc.py` plots the mean of all pairwise crossings, which differs by at most 1.16 K; both give the maximum of the cubic fit at ε = −2.87 %. The scripts were run in the original flat working directory on TRUBA, so some default paths point there; all of them take the input location as a command-line argument. Comments and printed labels are partly in Turkish (eksik = missing, sapma = deviation, DOGRULANDI = verified, artik = residual, tepe = peak).

## Requirements

Python ≥ 3.8 with NumPy (Matplotlib for the figures).

## Licence

Code: MIT (see `LICENSE`). Data: CC BY 4.0.

## Contact

Nejdet Paran, Department of Physics, Faculty of Science, Erciyes University, 38039 Kayseri, Türkiye — nejdetparan@erciyes.edu.tr — ORCID 0000-0003-0317-0438
