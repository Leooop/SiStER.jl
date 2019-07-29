#### [n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,quad,x,y,args)
#### interpolates marker properties to shear nodes
#### First cut - J.A. Olive, March 2011
## Modified by E. Mittelstaedt, April 2011, to allow variable inputs.
## Modified by B.Z. Klein, Spring 2014, for speedup
### Modified by B.Z. Klein, Summer 2014, for further speedup (vectorized)

# mutable struct DataInterp{T<:Real}
#     data::Array{T,2}
# end

function interp_markers_to_shear_nodes_parallel(xm,ym,icn,jcn,quad,x,y,args...)

    Nx=length(x)
    Ny=length(y)
    dx=diff(x)
    dy=diff(y)

    ### MITTELSTAEDT - check for number of properties to interpolate
    numV = length(args)


    ### MITTELSTAEDT ### establish interpolants matrices
    #n2interp = repmat(struct('data', zeros(Ny,Nx)), 1, numV);
    n2interp = Vector{Matrix{Float64}}(undef,numV)


    # Create interpolation object for the index with nearnest (Constant()) parametrization
    itpX = interpolate((x,), 1:length(x), Gridded(Constant()))
    itpY = interpolate((y,), 1:length(y), Gridded(Constant()))

    # Allow extrapolation with values outside the boundaries being equal to the boundaries values
    itpX = extrapolate(itpX, Flat())
    itpY = extrapolate(itpY, Flat())

    # Interpolate to {xm,ym}
    ## these are like the jcn and icn elsewhere, except the nodes are centered instead of upper left.
    ## this makes a lot of the indexing much simpler below.
    JCN = Int.(itpX(xm))
    ICN = Int.(itpY(ym))


    ###### Interior Cells

    center = (jcn.>1) .& (jcn.<Nx) .& (icn.>1) .& (icn.<Ny)
    shiftLeft = (jcn.<Nx-1) .& (icn.>1) .& (icn.<Ny)
    shiftUp = (jcn.>1) .& (jcn.<Nx) .& (icn.<Ny-1)
    shiftBoth = (jcn.<Nx-1) .& (icn.<Ny-1)

    cell1 = center .& ((xm.-x[JCN]) .> 0) .& ((ym .- y[ICN]) .> 0)  ## these are logical arrays that index the original quadrants
    cell2 = shiftLeft .& ((xm.-x[JCN]) .< 0) .& ((ym .- y[ICN]) .> 0)
    cell3 = shiftBoth .& ((xm.-x[JCN]) .< 0) .& ((ym .- y[ICN]) .< 0)
    cell4 = shiftUp .& ((xm.-x[JCN]) .> 0) .& ((ym .- y[ICN]) .< 0)

    ######### WEIGHTING (equal for now because that is what I'm running)

    wc1 = 0.25
    wc2 = 0.25
    wc3 = 0.25
    wc4 = 0.25


    ### cell 1 (i,j,1)

    dxm = xm[cell1] .- x[JCN[cell1]]
    dym = ym[cell1] .- y[ICN[cell1]]
    ddx = dx[JCN[cell1]]
    ddy = dy[ICN[cell1]]

    wm1 = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx.*ddy)
    #w1 = accumarray([ICN(cell1)', JCN(cell1)'], wm1, [Ny, Nx]) :
    w1 = accumarray(ICN[cell1], JCN[cell1], wm1, (Ny,Nx))
    # cell 2 (i, j-1, 2)

    dxm = xm[cell2] .- x[JCN[cell2]]
    dym = ym[cell2] .- y[ICN[cell2]]
    ddx = dx[JCN[cell2].-1]
    ddy = dy[ICN[cell2]]

    wm2 = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx.*ddy)
    #w2 = accumarray([ICN(cell2)', JCN(cell2)'], wm2, [Ny, Nx]);
    w2 = accumarray(ICN[cell2], JCN[cell2], wm2, (Ny,Nx))
    # cell 3 (i-1, j-1, 3)

    dxm = xm[cell3] .- x[JCN[cell3]]
    dym = ym[cell3] .- y[ICN[cell3]]
    ddx = dx[JCN[cell3].-1]
    ddy = dy[ICN[cell3].-1]

    wm3 = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx.*ddy)
    #w3 = accumarray([ICN(cell3)', JCN(cell3)'], wm3, [Ny, Nx])
    w3 = accumarray(ICN[cell3], JCN[cell3], wm3, (Ny,Nx))
    # cell 4 (i-1, j, 4)

    dxm = xm[cell4] .- x[JCN[cell4]]
    dym = ym[cell4] .- y[ICN[cell4]]
    ddx = dx[JCN[cell4]]
    ddy = dy[ICN[cell4].-1]

    wm4 = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx.*ddy)
    #w4 = accumarray([ICN(cell4)', JCN(cell4)'], wm4, [Ny, Nx]);
    w4 = accumarray(ICN[cell4], JCN[cell4], wm4, (Ny,Nx))
    ###loop over material properties to interpolate

    for vn = 1:numV

        w1_term = @spawn accumarray(ICN[cell1], JCN[cell1],
            args[vn][cell1].*wm1, (Ny,Nx))
        w2_term = @spawn accumarray(ICN[cell2], JCN[cell2],
            args[vn][cell2].*wm2, (Ny,Nx))
        w3_term = @spawn accumarray(ICN[cell3], JCN[cell3],
            args[vn][cell3].*wm3, (Ny,Nx))
        w4_term = @spawn accumarray(ICN[cell4], JCN[cell4],
            args[vn][cell4].*wm4, (Ny,Nx))



        n2interp[vn] = ((wc1.*fetch(w1_term))./w1 +
                            (wc2.*fetch(w2_term))./w2 +
                            (wc3.*fetch(w3_term))./w3 +
                            (wc4.*fetch(w4_term))./w4 )./ (wc1+wc2+wc4+wc4)

        # n2interp(vn).data = (wc1*accumarray([ICN(cell1)', JCN(cell1)'], args{vn}(cell1).*wm1)./w1 + ...
        #     wc2*accumarray([ICN(cell2)', JCN(cell2)'], args{vn}(cell2).*wm2)./w2 + ...
        #     wc3*accumarray([ICN(cell3)', JCN(cell3)'], args{vn}(cell3).*wm3)./w3 + ...
        #     wc4*accumarray([ICN(cell4)', JCN(cell4)'], args{vn}(cell4).*wm4)./w4)./...
        #     (wc1+wc2+wc4+wc4);
    end



    ###### EDGES

    ######### top edge

    topEdge = (jcn.>1) .& (jcn.<Nx) .& (icn.==1)
    shifted = (jcn.<Nx-1) .& (icn.==1)

    # cell 1

    cell1 = shifted .& (quad.==2)

    ddx = dx[JCN[cell1].-1]    # ???????
    ddy = dy[1]
    dxm = xm[cell1] .- x[JCN[cell1]]  # ????????
    dym = ym[cell1] .- y[ICN[cell1]]
    wm1 = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx.*ddy)
    #w1 = accumarray([ICN(cell1)', JCN(cell1)'], wm1, [1, Nx]);
    w1 = accumarray(ICN[cell1], JCN[cell1], wm1, (1,Nx))

    # cell 2

    cell2 = topEdge .& (quad.==1)

    ddx = dx[JCN[cell2]]
    ddy = dy[1]
    dxm = xm[cell2] .- x[JCN[cell2]]
    dym = ym[cell2] .- y[ICN[cell2]]
    wm2 = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx.*ddy)
    #w2  = accumarray([ICN(cell2)', JCN(cell2)'], wm2, [1, Nx]);
    w2 = accumarray(ICN[cell2], JCN[cell2], wm2, (1,Nx))

    #loop over material properties to interpolate

    for vn = 1:numV

        w1_term = @spawn accumarray(ICN[cell1], JCN[cell1], args[vn][cell1].*wm1, (1,Nx))
        w2_term = @spawn accumarray(ICN[cell2], JCN[cell2], args[vn][cell2].*wm2, (1,Nx))

        n2interp[vn][1,:] = ((wc1.*fetch(w1_term))./w1 +
                        (wc2.*fetch(w2_term))./w2) ./ (wc1+wc2)

        # temp = (wc1*accumarray([ICN(cell1)', JCN(cell1)'], args{vn}(cell1).*wm1)./w1 + ...
        #     wc2*accumarray([ICN(cell2)', JCN(cell2)'], args{vn}(cell2).*wm2)./w2)/...
        #     (wc1+wc2);
        # n2interp(vn).data(1,2:end) = temp(2:end);
    end

    ###### bottom edge

    bottomEdge = (jcn.>1) .& (jcn.<Nx) .& (icn.==(Ny-1))
    shifted = (jcn.<(Nx-1)) .& (icn.==(Ny-1))

    # cell 1

    cell1 = shifted .& (quad.==3)

    ddx = dx[JCN[cell1].-1]
    ddy = dy[Ny-1]
    dxm = xm[cell1] .- x[JCN[cell1]]
    dym = ym[cell1] .- y[end-1]
    wm1 = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx.*ddy)
    #w1 = accumarray([ones(sum(cell1),1), JCN(cell1)'], wm1, [1, Nx]);
    w1 = accumarray(ones(Int,sum(cell1)), JCN[cell1], wm1, (1,Nx))
    # cell 2

    cell2 = bottomEdge .& (quad.==4)

    ddx = dx[JCN[cell2]]
    ddy = dy[Ny-1]
    dxm = xm[cell2] .- x[JCN[cell2]]
    dym = ym[cell2] .- y[end-1]
    wm2 = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx.*ddy)
    #w2  = accumarray([ones(sum(cell2),1), JCN(cell2)'], wm2, [1, Nx]);
    w1 = accumarray(ones(Int,sum(cell2)), JCN[cell2], wm2, (1,Nx))

    ##loop over material properties to interpolate

    for vn = 1:numV

        w1_term = @spawn accumarray(ones(Int,sum(cell1)), JCN[cell1], args[vn][cell1].*wm1, (1,Nx))
        w2_term = @spawn accumarray(ones(Int,sum(cell2)), JCN[cell2], args[vn][cell2].*wm2, (1,Nx))

        n2interp[vn][Ny,:] = ((wc1.*fetch(w1_term))./w1 .+
                        (wc2.*fetch(w2_term))./w2) ./ (wc1+wc2)

        # temp = (wc1*accumarray([ones(sum(cell1),1), JCN(cell1)'], args{vn}(cell1).*wm1)./w1 + ...
        #     wc2*accumarray([ones(sum(cell2),1), JCN(cell2)'], args{vn}(cell2).*wm2)./w2)/...
        #     (wc1+wc2);
        # n2interp(vn).data(Ny,2:end) = temp(2:end);
    end

    ###### left edge

    leftEdge = (jcn.==1) .& (icn.>1) .& (icn.<Ny)
    shifted  = (jcn.==1) .& (icn.<(Ny-1))

    # cell 1

    cell1 = shifted .& (quad.==4)

    ddx = dx[1]
    ddy = dy[ICN[cell1].-1]
    dxm = xm[cell1] .- x[1]
    dym = ym[cell1] .- y[ICN[cell1]]
    wm1 = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx.*ddy)
    #w1 = accumarray([ICN(cell1)', ones(sum(cell1),1)], wm1, [Ny, 1]);
    w1 = accumarray(ICN[cell1], ones(Int,sum(cell1)), wm1, (Ny,1))

    # cell 2

    cell2 = leftEdge .& (quad.==1)

    ddx = dx[1]
    ddy = dy[ICN[cell2]]
    dxm = xm[cell2] .- x[1]
    dym = ym[cell2] .- y[ICN[cell2]]
    wm2 = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx.*ddy)
    #w2 = accumarray([ICN(cell2)', ones(sum(cell2),1)], wm2, [Ny, 1]);
    w2 = accumarray(ICN[cell2], ones(Int,sum(cell2)), wm2, (Ny,1))

    ##loop over material properties to interpolate

    for vn = 1:numV

        w1_term = @spawn accumarray(ICN[cell1], ones(Int,sum(cell1)), args[vn][cell1].*wm1, (Ny,1))
        w2_term = @spawn accumarray(ICN[cell2], ones(Int,sum(cell2)), args[vn][cell2].*wm2, (Ny,1))

        n2interp[vn][:, 1] = ((wc1.*fetch(w1_term))./w1 +
                        (wc2.*fetch(w2_term))./w2) ./ (wc1+wc2)

        # temp = (wc1*accumarray([ICN(cell1)', ones(sum(cell1),1)], args{vn}(cell1).*wm1)./w1 + ...
        #     wc2*accumarray([ICN(cell2)', ones(sum(cell2),1)], args{vn}(cell2).*wm2)./w2)/...
        #     (wc1+wc2);
        # n2interp(vn).data(2:end-1, 1) = temp(2:end);
    end

    ###### right edge

    rightEdge = (jcn.==Nx-1) .& (icn.>1) .& (icn.<Ny)
    shifted =   (jcn.==(Nx-1)) .& (icn.<(Ny-1))

    # cell 1

    cell1 = shifted .& (quad.==3)

    ddx = dx[Nx-1]
    ddy = dy[ICN[cell1].-1]
    dxm = xm[cell1] .- x[Nx-1]
    dym = ym[cell1] .- y[ICN[cell1]]
    wm1 = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx.*ddy)
    #w1 = accumarray([ICN(cell1)', ones(sum(cell1),1)], wm1, [Ny, 1]);
    w1 = accumarray(ICN[cell1], ones(Int,sum(cell1)), wm1, (Ny,1))
    # cell 2

    cell2 = rightEdge .& (quad.==2)

    ddx = dx[Nx-1]
    ddy = dy[ICN[cell2]]
    dxm = xm[cell2] .- x[Nx-1]
    dym = ym[cell2] .- y[ICN[cell2]]
    wm2 = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx.*ddy)
    #w2 = accumarray([ICN(cell2)', ones(sum(cell2),1)], wm2, [Ny, 1]);
    w2 = accumarray(ICN[cell2], ones(Int,sum(cell2)), wm2, (Ny,1))

    ##loop over material properties to interpolate

    for vn = 1:numV

        w1_term = @spawn accumarray(ICN[cell1], ones(Int,sum(cell1)), args[vn][cell1].*wm1, (Ny,1))
        w2_term = @spawn accumarray(ICN[cell2], ones(Int,sum(cell2)), args[vn][cell2].*wm2, (Ny,1))

        n2interp[vn][:,Nx] = ((wc1.*fetch(w1_term))./w1 +
                        (wc2.*fetch(w2_term))./w2) ./ (wc1+wc2)

        # temp = (wc1*accumarray([ICN(cell1)', ones(sum(cell1),1)], args{vn}(cell1).*wm1)./w1 + ...
        #     wc2*accumarray([ICN(cell2)', ones(sum(cell2),1)], args{vn}(cell2).*wm2)./w2)/...
        #     (wc1+wc2);
        # n2interp(vn).data(2:end-1, Nx) = temp(2:end);
    end

    #### CORNERS

    ## upper left

    upperLeft = (jcn.==1) .& (icn.==1) .& (quad.==1)

    ddx = dx[1]
    ddy = dy[1]
    dxm = xm[upperLeft] .- x[1]
    dym = ym[upperLeft] .- y[1]
    wm  = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx*ddy)
    wco = sum(wm)

    for vn = 1:numV
        n2interp[vn][1,1] = sum(args[vn][upperLeft].*wm)./wco
    end

    ## upper right

    upperRight = (icn.==1) .& (jcn.==(Nx-1)) .& (quad.==2)

    ddx = dx[Nx-1]
    ddy = dy[1]
    dxm = xm[upperRight] .- x[Nx-1]
    dym = ym[upperRight] .- y[1]
    wm  = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx*ddy)
    wco = sum(wm)

    for vn = 1:numV
        n2interp[vn][1,Nx] = sum(args[vn][upperRight].*wm)./wco
    end

    ## lower Right

    lowerRight = (icn.==Ny-1) .& (jcn.==(Nx-1)) .& (quad.==3)

    ddx = dx[Nx-1]
    ddy = dy[Ny-1]
    dxm = xm[lowerRight] .- x[Nx-1]
    dym = ym[lowerRight] .- y[Ny-1]
    wm  = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx*ddy)
    wco = sum(wm)

    for vn = 1:numV
        n2interp[vn][Ny,Nx] = sum(args[vn][lowerRight].*wm)./wco
    end

    # lower left

    lowerLeft = (icn.==(Ny-1)) .& (jcn.==1) .& (quad.==4)

    ddx = dx[1]
    ddy = dy[Ny-1]
    dxm = xm[lowerLeft] .- x[1]
    dym = ym[lowerLeft] .- y[Ny-1]
    wm  = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx*ddy)
    wco = sum(wm)

    for vn = 1:numV
        n2interp[vn][Ny,1] = sum(args[vn][lowerLeft].*wm)./wco
    end

    return [PropertyDict(Dict(:data => n2interp[i])) for i in eachindex(n2interp)]

end
