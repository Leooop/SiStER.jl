# initialize_marker_phases(Nphase,GEOM,xm,ym)
# this is where the identity of each marker (e.g., air, crust, mantle...)
# gets assigned following the geometry specified in the input file
# (flat layers, circle or rectangle)
# an alternative geometry can be added here.

function initialize_marker_phases(PARAMS,GEOM,xm,ym)

    # assign material identity on markers
    phm=Vector{Int}(undef,size(xm))

    for kk = 1:PARAMS.Nphase

        if GEOM[kk].type == 1 # layer

            phm[(ym.>=GEOM[kk].top) .& (ym.<GEOM[kk].bot)] .= kk

        elseif GEOM[kk].type == 2 # circular inclusion

            rm = sqrt.((xm.-GEOM[kk].x0).^2 .+ (ym.-GEOM[kk].y0).^2)
            phm[rm.<=GEOM[kk].rad] .= kk

        elseif GEOM[kk].type == 3 # rectangle

            phm[(ym.>=GEOM(kk).top) .& (ym.<GEOM[kk].bot) .& (xm.>=GEOM[kk].left) .& (xm.<GEOM[kk].right)] .= kk

        end

    end

    return phm

end

# Visualization :
# scatter(reshape(X,length(X)),reshape(Y,length(Y)))
# for i = 1:3
#     inds = (phm .== i)
#     scatter!(xm[inds],ym[inds],ms = 2)
# end
# scatter!()
