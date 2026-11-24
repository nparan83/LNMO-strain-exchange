import re
t = open('eps0_cfg0.in').read()
sp = []
for l in t.split('ATOMIC_SPECIES')[1].split('\n')[1:]:
    q = l.split()
    if len(q) >= 3 and re.match(r'^[A-Za-z]{1,2}\d*$', q[0]): sp.append(q[0])
    elif sp: break
rest = t.split('ATOMIC_POSITIONS')[1].split('\n',1)[1]
F = {}
for l in rest.split('\n'):
    q = l.split()
    if len(q) >= 4 and re.match(r'^[A-Za-z]{1,2}\d*$', q[0]): F[q[0]] = [float(x) for x in q[1:4]]
    elif F: break
M = {}
for L in F:
    if not L.startswith(('Mn','Ni')): continue
    f1, f2, f3 = F[L]
    m3 = int(round(2*f3))
    m1 = int(round(2*f1 - 0.5 - 0.5*m3)); m2 = int(round(2*f2 - 0.5 - 0.5*m3))
    M[L] = (m1, m2, m3)
    print("%-4s f=(%.3f,%.3f,%.3f)  m=(%+d,%+d,%+d)  parite=%d" % (L,f1,f2,f3,m1,m2,m3,sum((m1,m2,m3))%2))
pm = {sum(M[L])%2 for L in M if L.startswith('Mn')}
pn = {sum(M[L])%2 for L in M if L.startswith('Ni')}
print("rocksalt: Mn%s Ni%s -> %s" % (pm, pn, "OK" if len(pm)==1==len(pn) and pm!=pn else "BOZUK"))
def mk(name, f):
    out = t.replace("prefix='eps0_cfg0'", "prefix='eps0_%s'" % name); pat=''; Mt=0
    for k,s in enumerate(sp):
        if not s.startswith(('Mn','Ni')): continue
        g = f(*M[s]); mg = 0.70 if s.startswith('Mn') else 0.50
        pat += '+' if g>0 else '-'; Mt += g*(3 if s.startswith('Mn') else 2)
        out = re.sub(r'starting_magnetization\(\s*%d\s*\)\s*=\s*\S+' % (k+1),
                     'starting_magnetization(%d)=%s%.2f' % (k+1, '+' if g>0 else '-', mg), out)
    open('eps0_%s.in' % name,'w').write(out); print("%-5s %s  Mtot=%+d" % (name, pat, Mt))
mk('cfgG', lambda a,b,c: 1 if (a+b+c)%2==0 else -1)
mk('cfgA', lambda a,b,c: 1 if c%2==0 else -1)
mk('cfgC', lambda a,b,c: 1 if (a+b)%2==0 else -1)
