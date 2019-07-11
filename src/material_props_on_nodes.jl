#-------------------------------------------------------------------------
# Get material properties on nodes from advected properties of markers
# G.Ito 8/16 - replaces previous version, which relied on marker viscosity
# for Picard iterations (J.-A.O.)
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
n2interp = interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,rhom)
rho  = n2interp[1].data

# GET MARKER ELASTIC PROPERTIES  G.Ito 8/16
Gm=get_elastic_moduli(phm,MAT)
# pass shear modulus to nodes
n2interp = interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,1.0 ./Gm)
Gn = 1 ./(n2interp[1].data)
n2interp = interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,1.0 ./Gm);
Gs = 1 ./(n2interp[1].data)

# testn = interp_markers_to_normal_nodes(vec(vars["xm"]),vec(vars["ym"]),Int.(vec(vars["icn"])),Int.(vec(vars["jcn"])),vec(vars["x"]),vec(vars["y"]),1 ./vec(vars["Gm"]))
# Gnt=1 ./(testn[1].data)
#
# tests = interp_markers_to_shear_nodes(vec(vars["xm"]),vec(vars["ym"]),Int.(vec(vars["icn"])),Int.(vec(vars["jcn"])),Int.(vec(vars["qd"])),vec(vars["x"]),vec(vars["y"]),1 ./vec(vars["Gm"]))
# Gst = 1 ./(tests[1].data)

# PROPERTIES FOR PLASTICITY  G.Ito 8/16
cohes = get_cohesion(phm,ep,MAT) # cohesion depends on plastic strain
n2interp = interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,cohes)
Cohes_n = n2interp[1].data
n2interp = interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,cohes)
Cohes_s = n2interp[1].data

# GET FRICTION BASED ON MARKERS J.A. Olive 4/17
fric = get_friction(phm,ep,MAT) # friction depends on plastic strain
n2interp = interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,fric)
Mu_n = n2interp[1].data
n2interp = interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,fric)
Mu_s = n2interp[1].data

# ADVECTED strain rate invariant G.Ito 8/16
n2interp = interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,epsIIm)
epsII_n = n2interp[1].data
n2interp = interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,epsIIm)
epsII_s = n2interp[1].data


# OLD STRESSES AND PRESSURES G.Ito 8/16
n2interp = interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,sxxm)
#sxxOLD[:,:] = n2interp[1].data #Does not work :
sxxOLD = n2interp[1].data
sxxOLD_s = interp_normal_to_shear_nodes(sxxOLD,dx,dy)

n2interp = interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,sxym)
#sxyOLD(:,:) = n2interp(1).data;
sxyOLD = n2interp[1].data
sxyOLD_n = interp_shear_to_normal_nodes(sxyOLD)

#MIGHT WANT TO ADVECT PRESSURES (FOR SPEED?) G.Ito 8/16
pold = p
ps_old = interp_normal_to_shear_nodes(p,dx,dy)

EXYOLD = EXY
EXXOLD = EXX
EXX_sOLD = interp_normal_to_shear_nodes(EXX,dx,dy)
EXY_nOLD = interp_shear_to_normal_nodes(EXY)

#TEMPERATURE ARRAYS NEEDED FOR DUCTILE RHEOLOGY  G.Ito 8/16
n2interp =interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,Tm)
Ts = n2interp[1].data
n2interp = interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,Tm)
Tn=n2interp[1].data
