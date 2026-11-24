import numpy as np, re, sys
P = np.array([[1.,1.,0.],[1.,-1.,0.],[0.,0.,1.]]); Pi = np.linalg.inv(P)
TS = [np.zeros(3), np.array([1.,0.,0.])]
TPL = open('eps0_cfg0.in').read()
SPEC = [s for s in re.findall(r'^\s*([A-Za-z]{1,2}\d*)\s+\S+\s+\S+\s*$',
        TPL.split('ATOMIC_SPECIES')[1].split('ATOMIC_POSITIONS')[0], re.M)]
def rdin(p):
    t = open(p).read()
    cell = np.array([[float(x) for x in l.split()] for l in t.split('CELL_PARAMETERS')[1].split('\n')[1:4]])
    seg = t.split('ATOMIC_POSITIONS')[1]; head, rest = seg.split('\n',1)
    unit = head.strip().strip('()').strip(); lab=[]; fr=[]
    for l in rest.split('\n'):
        q = l.split()
        if len(q)>=4 and re.match(r'^[A-Za-z]{1,2}\d*$', q[0]):
            lab.append(re.sub(r'\d+$','',q[0])); fr.append([float(x) for x in q[1:4]])
        elif lab: break
    fr = np.array(fr)
    if not unit.startswith('cryst'): fr = fr @ np.linalg.inv(cell)
    return cell, lab, fr
def go(eps):
    cell, lab, fr = rdin('LNMO_Strain_Production/strain_%.2f/FM/scf_FM.in' % (eps/100.0))
    nc = P @ cell; L=[]; F=[]
    for t in TS:
        for s,f in zip(lab,fr): L.append(s); F.append(((f+t) @ Pi) % 1.0)
    F = np.array(F); cm=cn=0; tag=[]
    for s in L:
        if s=='Mn': cm+=1; tag.append('Mn%d'%cm)
        elif s=='Ni': cn+=1; tag.append('Ni%d'%cn)
        else: tag.append(s)
    M={}
    for g,f in zip(tag,F):
        if not g.startswith(('Mn','Ni')): continue
        n = [int(round(4*x)) % 4 for x in f]
        m3 = (n[2]//2) % 2
        for k in (0,1):
            if (n[k]-1-m3) % 2: sys.exit("eps=%d %s: yarim-tamsayi indeks, konum ideal degil: %s" % (eps,g,f))
        M[g] = ((n[0]-1-m3)//2, (n[1]-1-m3)//2, m3)
    pm = {sum(M[g])%2 for g in M if g.startswith('Mn')}; pn = {sum(M[g])%2 for g in M if g.startswith('Ni')}
    assert len(pm)==1==len(pn) and pm!=pn, "eps=%d rocksalt BOZUK Mn%s Ni%s" % (eps,pm,pn)
    pos = "ATOMIC_POSITIONS crystal\n" + "".join("%-4s %14.10f %14.10f %14.10f\n" % (g,*f) for g,f in zip(tag,F))
    cel = "CELL_PARAMETERS angstrom\n" + "".join("  %16.10f %16.10f %16.10f\n" % tuple(v) for v in nc)
    for nm, fn in [('FM', lambda a,b,c: 1), ('A', lambda a,b,c: 1 if c%2==0 else -1),
                   ('C', lambda a,b,c: 1 if (a+b)%2==0 else -1), ('G', lambda a,b,c: 1 if (a+b+c)%2==0 else -1)]:
        o = TPL
        o = re.sub(r'CELL_PARAMETERS[^\n]*\n(?:[^\n]*\n){3}', cel, o)
        o = re.sub(r'ATOMIC_POSITIONS[^\n]*\n(?:\s*[A-Za-z]{1,2}\d*\s+\S+\s+\S+\s+\S+[^\n]*\n)+', pos, o)
        o = o.replace("prefix='eps0_cfg0'", "prefix='e%+03d_%s'" % (eps,nm))
        pat=''
        for k,s in enumerate(SPEC):
            if not s.startswith(('Mn','Ni')): continue
            g = fn(*M[s]); mg = 0.70 if s.startswith('Mn') else 0.50; pat += '+' if g>0 else '-'
            o = re.sub(r'starting_magnetization\(\s*%d\s*\)\s*=\s*\S+' % (k+1),
                       'starting_magnetization(%d)=%s%.2f' % (k+1,'+' if g>0 else '-',mg), o)
        open('e%+03d_%s.in' % (eps,nm),'w').write(o)
        print("e%+03d_%s.in  %s  V=%.2f A^3  nat=%d" % (eps,nm,pat,abs(np.linalg.det(nc)),len(tag)))
for eps in (-4,-3,-2,-1,1,2,3,4): go(eps)
