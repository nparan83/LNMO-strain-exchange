import numpy as np, re
RY = 13605.693
MAG = [5,6,7,8,25,26,27,28]     # Mn1 Mn2 Ni1 Ni2 Mn3 Mn4 Ni3 Ni4 (atom sirasi)
EXP = {'FM':'++++++++', 'A':'+--++--+', 'C':'+-+-+-+-', 'G':'++--++--'}
def rd(f):
    t = open(f).read()
    ok = 'convergence has been achieved' in t
    e = float(re.findall(r'^!\s+total energy\s+=\s+(\S+)', t, re.M)[-1]) * RY
    b = t.split('Magnetic moment per site')[-1]; d = {}
    for l in b.split('\n'):
        g = re.match(r'\s+atom\s+(\d+).*magn=\s*(\S+)', l)
        if g: d[int(g.group(1))] = float(g.group(2))
    mt = re.findall(r'total magnetization\s+=\s+(\S+)', t)[-1]
    return ok, e, d, mt
files = [('FM','eps0_cfg0.out'), ('A','eps0_cfgA.out'), ('C','eps0_cfgC.out'), ('G','eps0_cfg1.out')]
M = np.array([[1.,-24.,-12.], [1.,-24.,12.], [1.,24.,-12.], [1.,24.,12.]])
Ev = []; bad = []
for l, f in files:
    ok, e, d, mt = rd(f); Ev.append(e)
    s = ''.join('+' if d.get(a,0) > 0 else '-' for a in MAG)
    good = ok and s == EXP[l]
    if not good: bad.append(l)
    print("%-3s conv=%s Mtot=%-6s desen=%s %s  |m|=%s" % (l, ok, mt, s, 'OK' if s==EXP[l] else 'BEKLENEN '+EXP[l],
          ' '.join('%+.2f' % d.get(a,0) for a in MAG)))
if bad: raise SystemExit("BOZUK: %s  -> fit yapilmadi" % bad)
Ev = np.array(Ev)
p, *_ = np.linalg.lstsq(M, Ev, rcond=None); res = Ev - M @ p
E0, Jp, Jq = p
print("\nJ1_par  = %+8.3f meV" % Jp)
print("J1_perp = %+8.3f meV" % Jq)
print("J||/J_| = %8.4f     [cos^2(theta) ongorusu: 1.010]" % (Jp/Jq))
print("<J1>    = %+8.3f meV     [izotropik pilot: 11.983]" % ((4*Jp+2*Jq)/6))
print("artik   = %s meV/hucre  (1 serbestlik derecesi)" % np.round(res,3))
print("\nBAGIMSIZ KONTRASTLAR:")
print("  J||  (E_C-E_FM)/48 = %+.3f   (E_G-E_A)/48 = %+.3f   fark = %.3f meV" %
      ((Ev[2]-Ev[0])/48, (Ev[3]-Ev[1])/48, abs((Ev[2]-Ev[0])/48-(Ev[3]-Ev[1])/48)))
print("  J_|  (E_A-E_FM)/24 = %+.3f   (E_G-E_C)/24 = %+.3f   fark = %.3f meV" %
      ((Ev[1]-Ev[0])/24, (Ev[3]-Ev[2])/24, abs((Ev[1]-Ev[0])/24-(Ev[3]-Ev[2])/24)))
