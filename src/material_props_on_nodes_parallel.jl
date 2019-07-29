#-------------------------------------------------------------------------
# Get material properties on nodes from advected properties of markers
# G.Ito 8/16 - replaces previous version, which relied on marker viscosity
# for Picard iterations (J.-A.O.)
# 20% slower than matlab version, mostly due to the absence of an equivalent to matlab's accumarray(). The implementation should be improved.
#-------------------------------------------------------------------------

# PHASE PROPORTIONS AT NORMAL AND SHEAR NODES. G.Ito 8/16
phase_n = interp_phases_to_normal_nodes(xm,ym,icn,jcn,x,y,phm,PARAMS)
phase_s = interp_phases_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,phm,PARAMS)
# phase_n and _s is a Ny*Nx*Nphase array containing the proportion
# of each phase at each node - this gets used in get_ductile_rheology
# functions

# OLD WAY TO INTERP PHASES: ONLY WORKED WELL WHEN MIXING 2 CONSECUTIVELY
# NUMBERED PHASES AT ANY NODE
# [n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,im);
# phase_n=n2interp(1).data;
# phase_n=round(phase_n*1e10)/1e10;  #prevents a case in which phase_n>NPhase
#
# [n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,im);
# phase_s=n2interp(1).data;
# phase_s=round(phase_s*1e10)/1e10; #prevents a case in which phase_n>NPhase

# GET MARKER DENSITIES
rhom = get_density(phm,Tm,MAT)
# pass density to nodes
rho = @spawn interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,rhom)[1].data



# GET MARKER ELASTIC PROPERTIES  G.Ito 8/16
Gm=get_elastic_moduli(phm,MAT)
# pass shear modulus to nodes
Gn = @spawn 1 ./(interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,1.0 ./Gm)[1].data)
Gs = @spawn 1 ./(interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,1.0 ./Gm)[1].data);

# testn = interp_markers_to_normal_nodes(vec(vars["xm"]),vec(vars["ym"]),Int.(vec(vars["icn"])),Int.(vec(vars["jcn"])),vec(vars["x"]),vec(vars["y"]),1 ./vec(vars["Gm"]))
# Gnt=1 ./(testn[1].data)
#
# tests = interp_markers_to_shear_nodes(vec(vars["xm"]),vec(vars["ym"]),Int.(vec(vars["icn"])),Int.(vec(vars["jcn"])),Int.(vec(vars["qd"])),vec(vars["x"]),vec(vars["y"]),1 ./vec(vars["Gm"]))
# Gst = 1 ./(tests[1].data)

# PROPERTIES FOR PLASTICITY  G.Ito 8/16
cohes = get_cohesion(phm,ep,MAT) # cohesion depends on plastic strain
Cohes_n = @spawn interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,cohes)[1].data
Cohes_s = @spawn interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,cohes)[1].data


# GET FRICTION BASED ON MARKERS J.A. Olive 4/17
fric = get_friction(phm,ep,MAT) # friction depends on plastic strain
Mu_n = @spawn interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,fric)[1].data
Mu_s = @spawn interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,fric)[1].data


# ADVECTED strain rate invariant G.Ito 8/16
epsII_n = @spawn interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,epsIIm)[1].data
epsII_s = @spawn interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,epsIIm)[1].data


# OLD STRESSES AND PRESSURES G.Ito 8/16
sxxOLD = @spawn interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,sxxm)[1].data
sxyOLD = @spawn interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,sxym)[1].data
sxxOLD = fetch(sxxOLD)
sxyOLD = fetch(sxyOLD)

sxxOLD_s = @spawn interp_normal_to_shear_nodes(sxxOLD,dx,dy)
sxyOLD_n = @spawn interp_shear_to_normal_nodes(sxyOLD)

#MIGHT WANT TO ADVECT PRESSURES (FOR SPEED?) G.Ito 8/16
pold = p
ps_old = @spawn interp_normal_to_shear_nodes(p,dx,dy)

EXYOLD = EXY
EXXOLD = EXX
EXX_sOLD = @spawn interp_normal_to_shear_nodes(EXX,dx,dy)
EXY_nOLD = @spawn interp_shear_to_normal_nodes(EXY)

#TEMPERATURE ARRAYS NEEDED FOR DUCTILE RHEOLOGY  G.Ito 8/16
Ts = @spawn interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,Tm)[1].data
Tn = @spawn interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,Tm)[1].data

rho = fetch(rho)
Gn = fetch(Gn)
Gs = fetch(Gs)
Cohes_n = fetch(Cohes_n)
Cohes_s = fetch(Cohes_s)
Mu_n = fetch(Mu_n)
Mu_s = fetch(Mu_s)
epsII_n = fetch(epsII_n)
epsII_s = fetch(epsII_s)
sxxOLD = fetch(sxxOLD)
sxxOLD_s = fetch(sxxOLD_s)
sxyOLD = fetch(sxyOLD)
sxyOLD_n = fetch(sxyOLD_n)
ps_old = fetch(ps_old)
EXX_sOLD = fetch(EXX_sOLD)
EXY_nOLD = fetch(EXY_nOLD)
Ts = fetch(Ts)
Tn = fetch(Tn)
