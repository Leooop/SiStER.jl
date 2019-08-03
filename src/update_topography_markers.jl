# update_topography_markers

# advect the marker chain that keeps track of topography
# in the current flow field
topo_x,topo_y = advect_markers(x,y,topo_x,topo_y,dx,dy,dt_m,vx,vy)

# locate the interface between sticky layer and left / right edge
if isempty(findfirst(x -> x.<0, topo_x))
    topoL=topo_y[1]
else
    #topoL=interp1(topo_x,topo_y,0)
    itpL = interpolate((topo_x,), topo_y, Gridded(Linear()))
    itpL = extrapolate(itpL, Line())
    topoL = itpL(0.0)
end

if isempty(findfirst(x -> x.>xsize, topo_x))
    topoR=topo_y[end]
else
    #topoR=interp1(topo_x,topo_y,xsize);
    itpR = interpolate((topo_x,), topo_y, Gridded(Linear()))
    itpR = extrapolate(itpR, Line())
    topoR = itpR(xsize)
end

# eliminate topography markers that left domain, keep the first one out on both sides
Iin=findall(0 .< topo_x .< xsize)

topo_x = topo_x[Iin]
pushfirst!(topo_x, 0.0)
push!(topo_x, xsize)

topo_y = topo_y[Iin]
pushfirst!(topo_y, topoL)
push!(topo_y, topoR)


if PARAMS.YNSurfaceProcesses==1
    # ERODE TOPOGRAPHY
    topo_y=topography_diffusion_solver(topo_x,topo_y,dt_m,PARAMS.topo_kappa)
    # RESET ROCK AND AIR (assumes topography is only interface between phase 1 and 2)
    #topomarkers=interp1(topo_x,topo_y,xm);
    itpM = interpolate((topo_x,), topo_y, Gridded(Linear()))
    topomarkers = itpM(xm)
    phm[(phm.==1) .& (ym.>=topomarkers)] .= 2
    phm[(phm.>=2) .& (ym.<topomarkers)] .= 1
end

# if there has been too much stretching, regrid the surface topography
if (maximum(diff(topo_x))>(5*topo_marker_spacing)) || (issorted(topo_x) == false)
    # surface regridding happens if somewhere 2 topo markers have been
    # stretched apart by more than 5 times the inital mean marker spacing
    # or if topo_x is no longer sorted due to compression.
    topo_xREGRID = range(0, xsize, length = Ntopo)
    #topo_yREGRID = interp1(topo_x,topo_y,topo_xREGRID(2:end-1));
    itpRG = interpolate((topo_x,), topo_y, Gridded(Linear()))
    itpRG = extrapolate(itpRG, Line())
    topo_yREGRID = itpRG(topo_xREGRID[2:end-1])
    topo_yREGRID=[topoL; topo_yREGRID; topoR]
    topo_x=collect(topo_xREGRID)
    topo_y=topo_yREGRID
    println("**REGRIDDING TOPOGRAPHY MARKERS**")
end
