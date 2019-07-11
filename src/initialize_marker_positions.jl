# SiStER_initialize_marker_positions(xsize,ysize,dx,dy,Mquad)
# assigns markers their inital position (coordinates)
# markers are seeded by cell quadrant, to make sure there is enough of them
# in each quadrant

function initialize_marker_positions(xsize,ysize,dx,dy,Mquad)

    # smallest quadrant sizes
    mdx=minimum(dx)/2
    mdy=minimum(dy)/2

    # creating a regular grid with step = that of the smallest quadrant
    xx=0:mdx:xsize
    yy=0:mdy:ysize

    nxx=length(xx)
    nyy=length(yy)


    midx=1
    xm = Float64[]
    ym = Float64[]
    for i=1:nyy-1
        for j=1:nxx-1

            xmtemp, ymtemp =seed_markers_uniformly(xx[j],yy[i],mdx,mdy,Mquad)
            # xm[midx:midx+Mquad-1]=xmtemp
            # ym[midx:midx+Mquad-1]=ymtemp
            append!(xm,xmtemp)
            append!(ym,ymtemp)
            midx=midx+Mquad

        end
    end

    return xm, ym
end

# visualization :

# inds = (phm .==1)
# scatter(X,Y)
# plt1 = scatter(X,Y,label="")
#     scatter!(xm[inds],ym[inds],ms = 0.5)
#     scatter!(xm[.~inds],ym[.~inds],ms = 0.5)
#     xlims!((0,5e3))
#     #ylims!((10e3-2.5e3,10e3+2.5e3))
#     ylims!((10e3-0.25e3,10e3+0.25e3))
