import numpy as np
kB = 0.0861733
def cross(A, B):
    T = A[:,0]; f = A[:,1]-B[:,1]; s = np.hypot(A[:,2], B[:,2])
    for i in range(len(f)-1):
        if f[i]*f[i+1] < 0:
            w = f[i]/(f[i]-f[i+1]); Tc = T[i] + w*(T[i+1]-T[i])
            sl = abs((f[i+1]-f[i])/(T[i+1]-T[i]))
            U = A[i,1] + w*(A[i+1,1]-A[i,1])
            return Tc, s[i:i+2].mean()/sl, U
    return None, None, None
J = np.loadtxt('J_strain.dat'); rows = []
print("eps   L-cifti kesisimleri (K)                          U4*     T_C(16/24)")
for eps in range(-5,6):
    try: d = np.loadtxt('mc_eps%d.dat' % eps)
    except Exception: print("%+3d  dosya yok" % eps); continue
    Ls = sorted(set(d[:,0].astype(int)))
    D = {}
    for L in Ls:
        x = d[d[:,0]==L][:,[1,3,4]]; D[L] = x[np.argsort(x[:,0])]
    txt = []; best = (None,None,None)
    for i in range(len(Ls)-1):
        for j in range(i+1, len(Ls)):
            Tc,e,U = cross(D[Ls[i]], D[Ls[j]])
            txt.append("%d/%d:%s" % (Ls[i], Ls[j], "%.1f" % Tc if Tc else "--"))
            if (Ls[i],Ls[j]) == (16,24): best = (Tc,e,U)
    r = J[np.argmin(np.abs(J[:,0]-eps))]
    Jp, Jq = r[1], r[2]
    Tmfa = (4*Jp+2*Jq)/3*1.5/kB
    print("%+3d  %-48s %s  %s" % (eps, " ".join(txt),
          "%.3f" % best[2] if best[2] else "  --  ",
          "%.1f +- %.1f" % (best[0], best[1]) if best[0] else "PENCERE DISI"))
    if best[0]: rows.append([eps, Jp, Jq, (4*Jp+2*Jq)/6, Tmfa, best[0], best[1]])
R = np.array(rows)
print("\n eps    J||      J_|    <J1>    T_MFA    T_C(MC)   T_C/T_C(0)  MC/MFA")
i0 = np.argmin(np.abs(R[:,0]))
for x in R:
    print("%+4d %8.3f %8.3f %7.3f %8.1f %9.1f %10.4f %8.4f"
          % (x[0], x[1], x[2], x[3], x[4], x[5], x[5]/R[i0,5], x[5]/x[4]))
np.savetxt('Tc_strain.dat', R, fmt='%.6f', header='eps J_par J_perp J_mean T_MFA T_C_MC err')
k = np.argmax(R[:,5])
print("\nT_C tepesi: eps = %+.0f%%  ->  %.1f K   |   T_C(0) = %.1f K   |   deney 273-280 K"
      % (R[k,0], R[k,5], R[i0,5]))
p = np.polyfit(R[:,0], R[:,5], 2)
print("parabol tepesi: eps = %+.2f%%   [<J1> tepesi -3.0%%, izotropi -0.74%%]" % (-p[1]/(2*p[0])))
