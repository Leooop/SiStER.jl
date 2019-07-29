# update elastic stresses on markers following a solve (but before advection)

ROT = get_rotation_rate(vx,vy,dx,dy,BC)
om=interp_shear_nodes_to_markers(ROT,x,y,xm,ym,icn,jcn)

# rotate markers
alpha = om.*dt_m
sxymtemp = sxxm.*sin.(2.0*alpha) .+ sxym.*cos.(2.0*alpha)
sxxm = sxxm.*(cos.(alpha).^2 .- sin.(alpha).^2) .-             sxym.*sin.(2.0*alpha)
sxym = sxymtemp
