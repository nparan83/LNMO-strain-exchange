import numpy as np, re
RY = 13605.693
MAG = [5,6,7,8,25,26,27,28]
EXP = {'FM':'++++++++','A':'+--++--+','C':'+-+-+-+-','G':'++--++--'}
M = np.array([[1.,-24.,-12.],[1.,-24.,12.],[1.,24.,-12.],[1.,24.,12.]])
def rd(f):
    t = open(f).read()
    ok = 'convergence has been achieved' in t
    e = float(re.findall(r'^!\s+total energy\s+=\s+(\S+)', t, re.M)[-1]) * RY
    b = t.split('Magnetic moment per site')[-1]; d = {}
    for l in b.split('\n'):
        g = re.match(r'\s+atom\s+(\d+).*magn=\s*(\S+)', l)
        if g: d[int(g.group(1))] = float(g.group(2))
    return ok, e, ''.join('+' if d.get(a,0)>0 else '-' for a in MAG), d
def files(eps):
    if eps == 0:
        return [('FM','eps0_cfg0.out'),('A','eps0_cfgA.out'),('C','eps0_cfgC.out'),('G','eps0_cfg1.out')]
    return [(l, 'e%+03d_%s.out' % (eps,l)) for l in ('FM','A','C','G')]
out = []
print("eps  durum        J||       J_|    J||/J_|     <J1>    artik    m(Mn)FM m(Ni)FM")
for eps in range(-5,6):
    E = []; prob = []
    for l,f in files(eps):
        try: ok,e,s,d = rd(f)
        except Exception as ex: prob.append("%s:%s" % (l,type(ex).__name__)); break
        if not ok: prob.append(l+":conv")
        if s != EXP[l]: prob.append("%s:desen(%s)" % (l,s))
        E.append(e)
        if l == 'FM': mMn, mNi = d[5], d[7]
    if len(E) < 4 or prob:
        print("%+3d  SORUN: %s" % (eps, ','.join(prob) if prob else 'dosya eksik')); continue
    E = np.array(E); p,*_ = np.linalg.lstsq(M,E,rcond=None); r = E - M@p
    Jp, Jq = p[1], p[2]; Jm = (4*Jp+2*Jq)/6
    print("%+3d  OK    %9.3f %9.3f %9.4f %9.3f %8.2f %8.2f %7.2f"
          % (eps, Jp, Jq, Jp/Jq, Jm, r[0], mMn, mNi))
    out.append([eps, Jp, Jq, Jm, r[0], mMn, mNi])
np.savetxt('J_strain.dat', np.array(out), fmt='%.6f',
           header='eps  J_par  J_perp  J_mean  res  m_Mn_FM  m_Ni_FM')
print("\nJ_strain.dat yazildi (%d nokta)" % len(out))
