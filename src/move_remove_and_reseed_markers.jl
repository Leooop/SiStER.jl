# SiStER Update Markers

xm_new,ym_new = advect_markers(x,y,xm,ym,dx,dy,dt_m,vx,vy)
xm=xm_new
ym=ym_new
# eliminate markers that left domain
Iin= findall((xm.<=xsize) .& (xm.>=0.0) .& (ym.>=0.0) .& (ym.<=ysize))

#ftest(xm,xsize) = @. xm<=xsize & xm>=0.0

#msg2='  markers removed: '
#msg=[msg2 num2str(length(xm)-length(Iin))]
#disp(msg)
xm=xm[Iin]
ym=ym[Iin]
phm=phm[Iin]
ep=ep[Iin]
epNH=epNH[Iin]
Tm=Tm[Iin]
idm=idm[Iin]
sxxm=sxxm[Iin]
sxym=sxym[Iin]
epsIIm=epsIIm[Iin]

# locate advected markers with respect to the eulerian grid
quad,icn,jcn = locate_markers_in_grid(xm,ym,x,y,dx,dy)

# check for holes in the marker distribution,
# patch with new markers if necessary
# those new markers immediately get assigned a value of phase (im), index
# (idm) and accumulated plastic strain (ep), i.e., the 2 variables that never get
# passed to nodes.
xm, ym, phm, Ifix, mp, ep, idm, Tm, sxxm, sxym, epNH, epsIIm = patch_marker_holes(icn,jcn,quad,Nx,Ny,Mquad,Mquad_crit,xm,ym,x,y,dx,dy,phm,ep,idm,Tm,sxxm,sxym,epNH, epsIIm)

# then they get assigned P, epsII and stresses from grid values

if minimum(Ifix)>0


    xmFIX = xm[Ifix]
    ymFIX = ym[Ifix]

    # pass temperature, pressure, strain rate and stresses to the new
    # markers from their nodal values
    # locate new markers with respect to the eulerian grid
    quadFIX,icnFIX,jcnFIX = locate_markers_in_grid(xmFIX,ymFIX,x,y,dx,dy)


    temp = interp_normal_nodes_to_markers(p,xc,yc,xmFIX,ymFIX,icnFIX,jcnFIX)
    pm[Ifix] = temp # pressure


end


# locate all markers with respect to the eulerian grid
qd,icn,jcn = locate_markers_in_grid(xm,ym,x,y,dx,dy)
