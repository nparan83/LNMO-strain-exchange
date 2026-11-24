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
    mt = re.findall(r'total magnetization\s+=\s+(\S+)', t)[-1]
    return ok, e, ''.join('+' if d.get(a,0)>0 else '-' for a in MAG), mt
rows = []
for eps, pre in [(-5,'e-05'), (0,'eps0'), (5,'e+05')]:
    fl = ([('FM','eps0_cfg0.out'),('A','eps0_cfgA.out'),('C','eps0_cfgC.out'),('G','eps0_cfg1.out')]
          if eps==0 else [(l,'%s_%s.out'%(pre,l)) for l in ('FM','A','C','G')])
    print("eps = %+d %%" % eps); E=[]; bad=0
    for l,f in fl:
        try: ok,e,s,mt = rd(f)
        except Exception as ex: print("   %-2s HATA %s" % (l,ex)); bad=1; break
        E.append(e); good = ok and s==EXP[l]; bad |= (not good)
        print("   %-2s %-12s conv=%s Mtot=%-6s %s %s" % (l,f,ok,mt,s,'OK' if s==EXP[l] else '<<< BEKLENEN '+EXP[l]))
    if bad: print("   -> atlandi\n"); continue
    E = np.array(E); p,*_ = np.linalg.lstsq(M,E,rcond=None); r = E - M@p
    rows.append((eps,p[1],p[2]))
    print("   J||=%.3f  J_|=%.3f  oran=%.4f  <J1>=%.3f  artik=%s\n"
          % (p[1],p[2],p[1]/p[2],(4*p[1]+2*p[2])/6, np.round(r,2)))
print(" eps      J||       J_|     J||/J_|     <J1>")
for e,x,y in rows: print("%+4d %9.3f %9.3f %9.4f %9.3f" % (e,x,y,x/y,(4*x+2*y)/6))
g = {-5:0.8640, 5:1.1285}
print("\nSENARYO TESTI  [eps=0: DFT 1.0538 vs geometri 1.0100]")
for e,x,y in rows:
    if e==0: continue
    print("  eps=%+d  DFT=%.4f  |  (a) sabit kayma: %.4f  |  (b) 5.4x buyutme: %.4f"
          % (e, x/y, g[e]+0.0438, 1+5.38*(g[e]-1)))
