# SiStER THERMAL SOLVE

# get previous temperature on nodes
if @isdefined(Told) == false
    Told = interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,Tm)[1].data
else
    Told[:,:] = interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,Tm)[1].data
end

# enforce Dirichlet boundary conditions to avoid mismatch between markers
# and nodes
if BCtherm.top[1]==1
    Told[1,:].=BCtherm.top[2]
end
if BCtherm.bot[1]==1
    Told[Ny,:].=BCtherm.bot[2]
end
if BCtherm.left[1]==1
    Told[:,1].=BCtherm.left[2]
end
if BCtherm.right[1]==1
    Told[:,Nx].=BCtherm.right[2]
end



# GET VARIABLE DIFFUSIVITY AND CP
if (haskey(MAT[1],:cp) || haskey(MAT[1],:k)) == false
    cpfield = fill(PARAMS.cpref,size(T))
    kfield = fill(PARAMS.kref,size(T))
    rhofield = fill(PARAMS.rhoref,size(T))
else
    km, cpm = get_thermal_properties(phm,MAT)
    rhom = get_density(phm,Tm,MAT)
    n2interp = interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,km,cpm,rhom)
    kfield=n2interp[1].data
    cpfield=n2interp[2].data
    rhofield=n2interp[3].data
end
# THERMAL SOLVE
T,_,_=thermal_solver_sparse_CFD(x,y,Told,rhofield,cpfield,kfield,dt_m,BCtherm,zeros(Float64,size(T)))



# temperature change
dT = T-Told
# enforce Dirichlet boundary conditions to avoid mismatch between markers
# and nodes
if BCtherm.top[1]==1
    dT[1,:] .=0.0
end
if BCtherm.bot[1]==1
    dT[Ny,:] .=0.0
end
if BCtherm.left[1]==1
    dT[:,1] .=0.0
end
if BCtherm.right[1]==1
    dT[:,Nx] .=0.0
end

Tm=interp_shear_nodes_to_markers(T,x,y,xm,ym,icn,jcn)

if PARAMS.ynTreset==1 # reset T=T0 in top layer
     Tm[phm.==1].=PARAMS.T0
end
