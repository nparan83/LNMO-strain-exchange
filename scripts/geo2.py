import numpy as np, re, itertools
IM = np.array(list(itertools.product([-1,0,1],repeat=3)))
def rd(p):
    t = open(p).read()
    b = t.split('Begin final coordinates')[-1].split('End final coordinates')[0] if 'Begin final coordinates' in t else t
    cl = b.split('CELL_PARAMETERS')[1].split('\n')[1:4]
    cell = np.array([[float(x) for x in l.split()] for l in cl])
    seg = b.split('ATOMIC_POSITIONS')[1]; head, rest = seg.split('\n',1)
    unit = head.strip().strip('()').strip(); sp=[]; fr=[]
    for l in rest.split('\n'):
        q = l.split()
        if len(q)>=4 and re.match(r'^[A-Za-z]{1,2}\d*$', q[0]):
            sp.append(re.sub(r'\d+$','',q[0])); fr.append([float(x) for x in q[1:4]])
        elif sp: break
    fr = np.array(fr); sp = np.array(sp)
    if not unit.startswith('cryst'): fr = fr @ np.linalg.inv(cell)
    return cell, sp, fr
def geom(cell, sp, fr):
    Ni = np.where(sp=='Ni')[0]; Mn = np.where(sp=='Mn')[0]; Ox = np.where(sp=='O')[0]
    B = []
    for i in Ni:
        for j in Mn:
            for s in IM:
                v = (fr[j]+s-fr[i]) @ cell; r = np.linalg.norm(v)
                if r > 4.4: continue
                best = None
                for k in Ox:
                    for so in IM:
                        v1 = (fr[k]+so-fr[i]) @ cell; r1 = np.linalg.norm(v1)
                        if r1 > 2.5: continue
                        v2 = (fr[k]+so-fr[j]-s) @ cell; r2 = np.linalg.norm(v2)
                        if r2 > 2.5: continue
                        th = np.degrees(np.arccos(np.clip((-v1)@(-v2)/(r1*r2),-1,1)))
                        if best is None or r1+r2 < best[0]: best = (r1+r2, r1, r2, th)
                if best: B.append((r, abs(v[2])/r, best[1], best[2], best[3]))
    return B
F = [(e, 'LNMO_Strain_Production/strain_%.2f/FM/scf_FM.in' % (e/100.0)) for e in range(-5,6)]
print("%4s %7s %7s %7s | %3s %3s | %8s %8s %7s | %8s %8s | %7s"
      % ("eps","a","b","V","n||","n_|","th_||","th_perp","dth","dNiO_||","dNiO_|","c2r"))
rows=[]
for e, f in F:
    try: cell, sp, fr = rd(f)
    except Exception as ex: print("%4d  HATA: %s (%s)" % (e, f, ex)); continue
    B = geom(cell, sp, fr)
    P = [b for b in B if b[1] < 0.5]; Q = [b for b in B if b[1] >= 0.5]
    if not P or not Q: print("%4d  siniflandirma basarisiz: n=%d" % (e,len(B))); continue
    tp = np.mean([b[4] for b in P]); tq = np.mean([b[4] for b in Q])
    dp = np.mean([b[2] for b in P]); dq = np.mean([b[2] for b in Q])
    V = abs(np.linalg.det(cell)); c2 = np.cos(np.radians(tp))**2/np.cos(np.radians(tq))**2
    print("%4d %7.4f %7.4f %7.2f | %3d %3d | %8.3f %8.3f %+7.3f | %8.4f %8.4f | %7.4f"
          % (e, np.linalg.norm(cell[0]), np.linalg.norm(cell[1]), V, len(P), len(Q), tp, tq, tp-tq, dp, dq, c2))
    rows.append([e,tp,tq,dp,dq,V])
np.savetxt('geo_theta.dat', np.array(rows), header="eps th_par th_perp dNiO_par dNiO_perp V", fmt="%.6f")
