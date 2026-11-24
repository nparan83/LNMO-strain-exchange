import numpy as np, sys
from multiprocessing import Pool
kB = 0.0861733
def run(arg):
    L, T, Jp, Jq, nEq, nMe, nBlk, seed = arg
    rng = np.random.default_rng(seed)
    ix, iy, iz = np.indices((L,L,L))
    par = (ix+iy+iz) % 2
    Smag = np.where(par==0, 1.5, 1.0)
    v = rng.normal(size=(L,L,L,3)); v /= np.linalg.norm(v,axis=-1,keepdims=True)
    s = v * Smag[...,None]
    Jax = (Jp, Jp, Jq); sub = [par==0, par==1]
    beta = 1.0/(kB*T); sig = 0.8
    def sweep(s, sig):
        na = nt = 0
        for m in sub:
            h = np.zeros_like(s)
            for ax,J in enumerate(Jax):
                h += J*(np.roll(s,1,axis=ax) + np.roll(s,-1,axis=ax))
            new = s + rng.normal(scale=sig, size=s.shape)*Smag[...,None]
            new *= (Smag/np.linalg.norm(new,axis=-1))[...,None]
            dE = -np.einsum('...i,...i->...', new-s, h)
            acc = (dE <= 0) | (rng.random(dE.shape) < np.exp(-beta*np.clip(dE,0,600/beta)))
            s = np.where((m & acc)[...,None], new, s)
            na += (m & acc).sum(); nt += m.sum()
        return s, na/nt
    for i in range(nEq):
        s, r = sweep(s, sig)
        if i % 200 == 199: sig = float(np.clip(sig*(1.6 if r > 0.5 else 0.7), 0.05, 3.0))
    Stot = Smag.sum(); nb = nMe//nBlk; U = []; Mb = []
    for _ in range(nBlk):
        m2 = m4 = mv = 0.0
        for _ in range(nb):
            s, _r = sweep(s, sig)
            mm = np.linalg.norm(s.sum(axis=(0,1,2)))/Stot
            mv += mm; m2 += mm*mm; m4 += mm**4
        m2 /= nb; m4 /= nb; mv /= nb
        U.append(1.0 - m4/(3*m2*m2)); Mb.append(mv)
    U = np.array(U); Mb = np.array(Mb)
    return L, T, Mb.mean(), U.mean(), U.std(ddof=1)/np.sqrt(nBlk)
if __name__ == '__main__':
    eps = int(sys.argv[1])
    D = np.loadtxt('J_strain.dat')
    Jp, Jq = D[np.argmin(np.abs(D[:,0]-eps))][1:3]
    Tmfa = (4*Jp+2*Jq)/3*1.5/kB
    Tc0 = 0.72*Tmfa
    Ts = np.round(np.linspace(Tc0-32, Tc0+32, 17), 1)
    jobs = [(L,T,Jp,Jq,20000,120000,12,10007*L+int(T*10)+eps) for L in (8,12,16,24) for T in Ts]
    print("# eps=%+d  J||=%.4f  J_|=%.4f meV  T_MFA=%.1f K  pencere %.1f-%.1f K"
          % (eps,Jp,Jq,Tmfa,Ts[0],Ts[-1]))
    print("# L        T        <m>        U4       err")
    with Pool(28) as p:
        for L,T,m,u,e in sorted(p.imap_unordered(run, jobs)):
            print("%4d %8.1f %10.5f %9.5f %9.5f" % (L,T,m,u,e), flush=True)
