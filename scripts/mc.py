import numpy as np, sys
from multiprocessing import Pool
kB = 0.0861733
def run(arg):
    L, T, Jp, Jq, nEq, nMe, seed = arg
    rng = np.random.default_rng(seed)
    ix, iy, iz = np.indices((L,L,L))
    par = (ix+iy+iz) % 2                 # 0 -> Mn (S=1.5), 1 -> Ni (S=1.0)
    Smag = np.where(par==0, 1.5, 1.0)
    v = rng.normal(size=(L,L,L,3)); v /= np.linalg.norm(v,axis=-1,keepdims=True)
    s = v * Smag[...,None]
    Jax = (Jp, Jp, Jq)
    sub = [par==0, par==1]
    beta = 1.0/(kB*T); sig = 0.8
    def field(s):
        h = np.zeros_like(s)
        for ax,J in enumerate(Jax):
            h += J*(np.roll(s,1,axis=ax) + np.roll(s,-1,axis=ax))
        return h
    def sweep(s):
        for m in sub:
            h = field(s)
            d = rng.normal(scale=sig, size=s.shape)
            new = s + d*Smag[...,None]
            new *= (Smag/np.linalg.norm(new,axis=-1))[...,None]
            dE = -np.einsum('...i,...i->...', new-s, h)
            acc = (dE <= 0) | (rng.random(dE.shape) < np.exp(-beta*np.clip(dE,0,700/beta)))
            upd = m & acc
            s = np.where(upd[...,None], new, s)
        return s
    for _ in range(nEq): s = sweep(s)
    Stot = Smag.sum(); m2 = m4 = 0.0
    for _ in range(nMe):
        s = sweep(s)
        mm = np.linalg.norm(s.sum(axis=(0,1,2)))/Stot
        m2 += mm*mm; m4 += mm**4
    m2 /= nMe; m4 /= nMe
    return L, T, np.sqrt(m2), 1.0 - m4/(3*m2*m2)
if __name__ == '__main__':
    eps = int(sys.argv[1])
    D = np.loadtxt('J_strain.dat')
    row = D[np.argmin(np.abs(D[:,0]-eps))]
    Jp, Jq = row[1], row[2]
    Ts = np.arange(200., 401., 10.)
    jobs = [(L,T,Jp,Jq,3000,8000,1000+7*L+int(T)) for L in (8,12,16) for T in Ts]
    print("# eps=%+d  J||=%.3f  J_|=%.3f meV" % (eps,Jp,Jq))
    print("# L      T        <m>       U4")
    with Pool(28) as p:
        for L,T,m,u in sorted(p.imap_unordered(run, jobs)):
            print("%4d %8.1f %10.5f %10.5f" % (L,T,m,u), flush=True)
