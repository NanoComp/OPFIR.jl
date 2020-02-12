using OPFIR

p = OPFIR.N2O(L=15., radius=0.25, pressure=50, power=0.25, num_layers=6, JL=14, n_vib=10,
       σ_DD=35, σ_GKC=15, σ_SPT=15,
       pbfactor=1., qclbroadening=1*1e6)

a, solnew, pnew, taus = OPFIR.outputpower(p, "U", "TE01", mumps_solver=0)
println(a[1])
