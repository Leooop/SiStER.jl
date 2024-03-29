# SiStER Initialize

#PARAMS.Nphase = Nphase; # for convenience ==> already added in input file

# construct staggered grids
X,Y,x,y,xc,yc,dx,dy,Nx,Ny = initialize_grid(xsize,ysize,GRID)

# initialize marker arrays and positions
xm, ym = initialize_marker_positions(xsize,ysize,dx,dy,Mquad)

# locate markers with respect to grid
qd,icn,jcn = locate_markers_in_grid(xm,ym,x,y,dx,dy)

# assign marker phases
phm = initialize_marker_phases(PARAMS,GEOM,xm,ym)

# initialize marker plastic strain (to zero) and strain rate (to one)
ep=zeros(Float64,size(xm))
epNH=copy(ep)
epsIIm=ones(Float64,size(xm))

# initialize marker stresses
sxxm=zeros(Float64,size(xm))
sxym=copy(sxxm)

# initialize marker index (a unique number to identify and track each marker)
idm=collect(1:length(xm))

# initialize temperature structure on nodes
T = PARAMS.a0 .+ PARAMS.a1.*Y .+ PARAMS.a2.*Y.^2 .+ PARAMS.a3.*Y.^3
T = T .+ PARAMS.amp.*sin.(2 .*pi.*X./PARAMS.lam)
if PARAMS.ynTreset==1 # reset T=T0 in top layer
    T[T.<PARAMS.T0].=PARAMS.T0
end
# pass initial nodal T to markers
Tm=interp_shear_nodes_to_markers(T,x,y,xm,ym,icn,jcn)
Tm0=copy(Tm)

# test :
# ext = extrema(Tm)
# inds = 1e-8 .< Tm .< 2e-8
# scatter(xm[inds],ym[inds],aspect_ratio = :equal)
# xlims!((6.5,9))




# initialize nodal strain rate and other useful arrays
EXX=zeros(Float64,size(X))
EXY=zeros(Float64,size(X))
vx=zeros(Float64,size(X))
vy=zeros(Float64,size(X))
v=copy(vx)
p=1e12.*ones(Float64,size(EXX))  #initialize to be high so plasticity doesnt activate at t=1, pit=1;
etan_new=zeros(Float64,Ny,Nx)

#-------------------------------------------------------------------------
# initialize dt_m small to keep things elastic & no plasticity at t=1, G.Ito
#-------------------------------------------------------------------------
if @isdefined(dt_m) == false
    dt_m = 1.0
end

# initialize marker chain to track base of layer 1 (sticky layer)
Ntopo=PARAMS.Ntopo_markers
topo_x=collect(range(0,xsize,length = Ntopo))
topo_y=GEOM[1].bot.*ones(size(topo_x))
topo_marker_spacing=mean(diff(topo_x)) # initial mean spacing of topography markers
