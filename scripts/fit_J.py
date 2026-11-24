import numpy as np, re, os
RY = 13605.693
d = np.load('design.npz', allow_pickle=True)
M, cfg, ms = d['M'], d['configs'], d['mag_sym']
MAG = [5,6,7,8,25,26,27,28]
E, bad = [], []
print("=== KONFIG DENETIMI ===")
for i in range(6):
    f = 'eps0_cfg%d.out' % i
    if not os.path.exists(f):
        print("cfg%d: DOSYA YOK" % i); bad.append(i); E.append(np.nan); continue
    t = open(f).read()
    conv = t.count('convergence has been achieved')
    e = re.findall(r'^!\s+total energy\s+=\s+(\S+)', t, re.M)
    mt = re.findall(r'total magnetization\s+=\s+(\S+)', t)
    blk = t.split('Magnetic moment per site')[-1].split('\n')
    mm = {}
    for ln in blk:
        g = re.match(r'\s+atom\s+(\d+).*magn=\s*(\S+)', ln)
        if g: mm[int(g.group(1))] = float(g.group(2))
    sg = [1 if mm.get(a,0) > 0 else -1 for a in MAG]
    ok = (sg == list(cfg[i]))
    if not ok or conv == 0: bad.append(i)
    E.append(float(e[-1]) * RY if e else np.nan)
    print("cfg%d conv=%d  E=%.6f Ry  Mtot=%s  desen=%s  |m|=%s" %
          (i, conv, float(e[-1]) if e else 0, mt[-1] if mt else '?',
           'OK' if ok else 'BOZUK ' + ''.join('+' if s>0 else '-' for s in sg),
           ' '.join('%.2f' % mm.get(a,0) for a in MAG)))
E = np.array(E)
print("\nBOZUK/EKSIK:", bad if bad else "yok")
if bad: raise SystemExit("Fit yapilmadi.")
p, *_ = np.linalg.lstsq(M, E, rcond=None)
res = E - M @ p
E0, J1, J2, J3 = p
print("\n=== FIT (H = -sum_<ij> J_ij S_i.S_j, J>0 = FM) ===")
for n, v in zip(['J1 Mn-Ni', 'J2 Ni-Ni', 'J3 Mn-Mn'], [J1, J2, J3]):
    print("  %s = %+8.3f meV = %+8.2f K" % (n, v, v / 0.0861733))
print("  E0 = %.4f meV" % E0)
print("\nartiklar (meV/hucre):", np.round(res, 3))
print("max|artik| = %.3f meV/hucre = %.3f meV/f.u." % (np.abs(res).max(), np.abs(res).max()/4))
print("|J2/J1| = %.3f   |J3/J1| = %.3f" % (abs(J2/J1), abs(J3/J1)))
print("\n=== CAPRAZ KONTROL: saf FiM (cfg1-cfg0) ===")
dJ1 = (E[1] - E[0]) / 72.0
print("  dE = %.2f meV/hucre = %.2f meV/f.u." % (E[1]-E[0], (E[1]-E[0])/4))
print("  J1(FiM tek basina) = %.3f meV   |   J1(fit) = %.3f meV" % (dJ1, J1))
print("  fark = %.1f %%" % (100*abs(dJ1-J1)/abs(J1)))
print("  paper-1: 192.2 -> 10.68 meV ; 215.96 -> 12.00 meV")
SA, SB, zAB, zAA, zBB = 1.0, 1.5, 6, 12, 12
A = (2/3)*np.array([[SA*(SA+1)*zAA*J2, np.sqrt(SA*(SA+1)*SB*(SB+1))*zAB*J1],
                    [np.sqrt(SA*(SA+1)*SB*(SB+1))*zAB*J1, SB*(SB+1)*zBB*J3]])
print("\nMFA T_C (ust sinir, referans) = %.0f K   [deney ~280 K]" % (np.linalg.eigvalsh(A).max()/0.0861733))
np.savez('J_eps0.npz', p=p, res=res, E=E)
