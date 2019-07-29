#-----------------------------------------------------------------------------
# Updates stresses on markers for CURRENT solution.  Stress rotation
# occurs after solutions are output
# G.Ito 8/16
#------------------------------------------------------------------------------

# Compute STRESS Changes on nodes, interpolate to markers, and apply to marker stresses
dsxx = (2.0.*etan.*EXX.-sxxOLD).*Zn
dsxy = (2.0.*etas.*EXY.-sxyOLD).*Zs

dsxxm = interp_normal_nodes_to_markers(dsxx,xc,yc,xm,ym,icn,jcn)
dsxym = interp_shear_nodes_to_markers(dsxy,x,y,xm,ym,icn,jcn)
sxxm = sxxm.+dsxxm
sxym = sxym.+dsxym
