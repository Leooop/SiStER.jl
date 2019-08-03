#   advect_markers(x,y,xm,ym,dx,dy,tstep,vx,vy)
#
# Uses a Fourth-Order Runge-Kutta method to advect markers
# in the current velocity field, over the advection time step tstep = dt_m
#
# First cut J.-A. Olive, 2011-2012 - modified by B.Z. Klein 2013-2014

function advect_markers(x,y,xm,ym,dx,dy,tstep,vx,vy)



    # point A (particle pts)
    XA = xm
    YA = ym
    qd,icn,jcn = locate_markers_in_grid(XA,YA,x,y,dx,dy)
    VxA,VyA = get_marker_velocities(qd,icn,jcn,x,y,xm,ym,vx,vy,dx,dy)

    # point B
    XB = @. XA + 0.5*tstep*VxA
    YB = @. YA + 0.5*tstep*VyA
    qd,icn,jcn = locate_markers_in_grid(XB,YB,x,y,dx,dy)
    VxB,VyB = get_marker_velocities(qd,icn,jcn,x,y,xm,ym,vx,vy,dx,dy)

    # point C
    XC = @. XA + 0.5*tstep*VxB
    YC = @. YA + 0.5*tstep*VyB
    qd,icn,jcn = locate_markers_in_grid(XC,YC,x,y,dx,dy)
    VxC,VyC = get_marker_velocities(qd,icn,jcn,x,y,xm,ym,vx,vy,dx,dy)

    # point D
    XD = @. XA + tstep*VxC
    YD = @. YA + tstep*VyC
    qd,icn,jcn = locate_markers_in_grid(XD,YD,x,y,dx,dy)
    VxD,VyD = get_marker_velocities(qd,icn,jcn,x,y,xm,ym,vx,vy,dx,dy)


    # effective velocity
    vx_eff = @. (1/6)*(VxA + 2*VxB + 2*VxC + VxD)
    vy_eff = @. (1/6)*(VyA + 2*VyB + 2*VyC + VyD)


    ## Calculate new coordinates
    xm_new = @. XA + tstep*vx_eff
    ym_new = @. YA + tstep*vy_eff

    return xm_new, ym_new, vx_eff, vy_eff

end
