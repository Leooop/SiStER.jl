# interp_shear_nodes_to_markers(varS,x,y,xm,ym,icn,jcn)
# interpolates properties from shear nodes to markers

function interp_shear_nodes_to_markers(varS,x,y,xm,ym,icn,jcn)

    m,n = size(varS)

    xnodes = [x[jcn] x[jcn.+1] x[jcn.+1] x[jcn]]
    ynodes = [y[icn] y[icn] y[icn.+1] y[icn.+1]]

    lin = LinearIndices(size(varS))
    INDEX = Vector{Int}(undef,length(icn))
    for i = 1:length(icn)
        INDEX[i] = lin[icn[i],jcn[i]]
    end

    VARnodes = [varS[INDEX] varS[INDEX.+m] varS[INDEX.+m.+1] varS[INDEX.+1]]

    varm = interp_grid_to_marker_vector(xnodes,ynodes,VARnodes,xm,ym)

    return varm
end
