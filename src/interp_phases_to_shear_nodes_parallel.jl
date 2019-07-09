# interp_phases_to_shear_nodes(xm,ym,icn,jcn,quad,x,y,phases,PARAMS.Nphase)
# B.Z. Klein, July 2017, an interp function specific to phase (to shear nodes), to enable
# exact mixing of several phases

function interp_phases_to_shear_nodes_parallel(xm,ym,icn,jcn,quad,x,y,phases,PARAMS)

    Nx=length(x)
    Ny=length(y)
    dx=diff(x)
    dy=diff(y)

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

    phasesS = zeros(Float64,(Ny, Nx, PARAMS.Nphase))

    ## Interior Cells

    center = (jcn.>1) .& (jcn.<Nx) .& (icn.>1) .& (icn.<Ny)
    shiftLeft = (jcn.<Nx-1) .& (icn.>1) .& (icn.<Ny)
    shiftUp = (jcn.>1) .& (jcn.<Nx) .& (icn.<Ny-1)
    shiftBoth = (jcn.<Nx-1) .& (icn.<Ny-1)


    cell1 = center .& ((xm.-x[JCN]) .> 0) .& ((ym .- y[ICN]) .> 0)  ## these are logical arrays that index the original quadrants
    cell2 = shiftLeft .& ((xm.-x[JCN]) .< 0) .& ((ym .- y[ICN]) .> 0)
    cell3 = shiftBoth .& ((xm.-x[JCN]) .< 0) .& ((ym .- y[ICN]) .< 0)
    cell4 = shiftUp .& ((xm.-x[JCN]) .> 0) .& ((ym .- y[ICN]) .< 0)

    ### WEIGHTING (equal for now because that is what I'm running)

    wc1 = 0.25
    wc2 = 0.25
    wc3 = 0.25
    wc4 = 0.25


    # cell 1 (i,j,1)

    dxm = xm[cell1] .- x[JCN[cell1]]
    dym = ym[cell1] .- y[ICN[cell1]]
    ddx = dx[JCN[cell1]]
    ddy = dy[ICN[cell1]]

    wm1 = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx.*ddy)
    #w1 = accumarray([ICN(cell1)', JCN(cell1)'], wm1, [Ny, Nx]) :
    w1 = zeros(Float64, (Ny,Nx))
    ICNp, JCNp = ICN[cell1], JCN[cell1]
    for i in eachindex(ICNp) # faster without multithreading
        w1[ICNp[i],JCNp[i]] += wm1[i]
    end
    # cell 2 (i, j-1, 2)

    dxm = xm[cell2] .- x[JCN[cell2]]
    dym = ym[cell2] .- y[ICN[cell2]]
    ddx = dx[JCN[cell2].-1]
    ddy = dy[ICN[cell2]]

    wm2 = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx.*ddy)
    #w2 = accumarray([ICN(cell2)', JCN(cell2)'], wm2, [Ny, Nx]);
    w2 = zeros(Float64, (Ny,Nx))
    ICNp, JCNp = ICN[cell2], JCN[cell2]
    for i in eachindex(ICNp)
        w2[ICNp[i],JCNp[i]] += wm2[i]
    end

    # cell 3 (i-1, j-1, 3)

    dxm = xm[cell3] .- x[JCN[cell3]]
    dym = ym[cell3] .- y[ICN[cell3]]
    ddx = dx[JCN[cell3].-1]
    ddy = dy[ICN[cell3].-1]

    wm3 = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx.*ddy)
    #w3 = accumarray([ICN(cell3)', JCN(cell3)'], wm3, [Ny, Nx])
    w3 = zeros(Float64, (Ny,Nx))
    ICNp, JCNp = ICN[cell3], JCN[cell3]
    for i in eachindex(ICNp)
        w3[ICNp[i],JCNp[i]] += wm3[i]
    end

    # cell 4 (i-1, j, 4)

    dxm = xm[cell4] .- x[JCN[cell4]]
    dym = ym[cell4] .- y[ICN[cell4]]
    ddx = dx[JCN[cell4]]
    ddy = dy[ICN[cell4].-1]

    wm4 = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx.*ddy)
    #w4 = accumarray([ICN(cell4)', JCN(cell4)'], wm4, [Ny, Nx]);
    w4 = zeros(Float64, (Ny,Nx))
    ICNp, JCNp = ICN[cell4], JCN[cell4]
    for i in eachindex(ICNp)
        w4[ICNp[i],JCNp[i]] += wm4[i]
    end

    #loop over material properties to interpolate

    Threads.@threads for n = 1:PARAMS.Nphase
        w1_term = zeros(Float64, (Ny,Nx))
        w2_term = copy(w1_term)
        w3_term = copy(w1_term)
        w4_term = copy(w1_term)

        phaseMask = (phases .== n)

        ICNp, JCNp, wm1p = ICN[cell1 .& phaseMask], JCN[cell1 .& phaseMask], wm1[phaseMask[cell1]]
        for i in eachindex(ICNp)
            w1_term[ICNp[i],JCNp[i]] += wm1p[i]
        end
        ICNp, JCNp, wm2p = ICN[cell2 .& phaseMask], JCN[cell2 .& phaseMask], wm2[phaseMask[cell2]]
        for i in eachindex(ICNp)
            w2_term[ICNp[i],JCNp[i]] += wm2p[i]
        end
        ICNp, JCNp, wm3p = ICN[cell3 .& phaseMask], JCN[cell3 .& phaseMask], wm3[phaseMask[cell3]]
        for i in eachindex(ICNp)
            w3_term[ICNp[i],JCNp[i]] += wm3p[i]
        end
        ICNp, JCNp, wm4p = ICN[cell4 .& phaseMask], JCN[cell4 .& phaseMask], wm4[phaseMask[cell4]]
        for i in eachindex(ICNp)
            w4_term[ICNp[i],JCNp[i]] += wm4p[i]
        end

        phasesS[:,:,n] = ((wc1.*w1_term)./w1 +
                        (wc2.*w2_term)./w2 +
                        (wc3.*w3_term)./w3 +
                        (wc4.*w4_term)./w4 )./ (wc1+wc2+wc4+wc4)

    end

    ## EDGES

    ### top edge

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
    w1 = zeros(Float64, (1,Nx))
    ICNp, JCNp = ICN[cell1], JCN[cell1]
    for i in eachindex(ICNp)
        w1[ICNp[i],JCNp[i]] += wm1[i]
    end
    # cell 2

    cell2 = topEdge .& (quad.==1)

    ddx = dx[JCN[cell2]]
    ddy = dy[1]
    dxm = xm[cell2] .- x[JCN[cell2]]
    dym = ym[cell2] .- y[ICN[cell2]]
    wm2 = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx.*ddy)
    #w2  = accumarray([ICN(cell2)', JCN(cell2)'], wm2, [1, Nx]);
    w2 = zeros(Float64, (1,Nx))
    ICNp, JCNp = ICN[cell2], JCN[cell2]
    for i in eachindex(ICNp)
        w2[ICNp[i],JCNp[i]] += wm2[i]
    end

    #loop over material properties to interpolate

    Threads.@threads for n = 1:PARAMS.Nphase

        w1_term = zeros(Float64, (1,Nx))
        w2_term = copy(w1_term)

        phaseMask = (phases .== n)

        ICNp, JCNp, wm1p = ICN[cell1 .& phaseMask], JCN[cell1 .& phaseMask], wm1[phaseMask[cell1]]
        for i in eachindex(ICNp)
            w1_term[ICNp[i],JCNp[i]] += wm1p[i]
        end
        ICNp, JCNp, wm2p = ICN[cell2 .& phaseMask], JCN[cell2 .& phaseMask], wm2[phaseMask[cell2]]
        for i in eachindex(ICNp)
            w2_term[ICNp[i],JCNp[i]] += wm2p[i]
        end

        phasesS[1,:,n] = ((wc1.*w1_term)./w1 +
                        (wc2.*w2_term)./w2) ./ (wc1+wc2)

    end

    ### bottom edge

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
    w1 = zeros(Float64, (1,Nx))
    ICNp, JCNp = ones(Int,sum(cell1)), JCN[cell1]
    for i in eachindex(ICNp)
        w1[ICNp[i],JCNp[i]] += wm1[i]
    end

    # cell 2

    cell2 = bottomEdge .& (quad.==4)

    ddx = dx[JCN[cell2]]
    ddy = dy[Ny-1]
    dxm = xm[cell2] .- x[JCN[cell2]]
    dym = ym[cell2] .- y[end-1]
    wm2 = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx.*ddy)
    #w2  = accumarray([ones(sum(cell2),1), JCN(cell2)'], wm2, [1, Nx]);
    w2 = zeros(Float64, (1,Nx))
    ICNp, JCNp = ones(Int,sum(cell2)), JCN[cell2]
    for i in eachindex(ICNp)
        w2[ICNp[i],JCNp[i]] += wm2[i]
    end
    #loop over material properties to interpolate

    Threads.@threads for n = 1:PARAMS.Nphase

        w1_term = zeros(Float64, (1,Nx))
        w2_term = copy(w1_term)

        phaseMask = (phases .== n)

        ICNp, JCNp, wm1p = ones(Int,sum(cell1 .& phaseMask)), JCN[cell1 .& phaseMask], wm1[phaseMask[cell1]]
        for i in eachindex(ICNp)
            w1_term[ICNp[i],JCNp[i]] += wm1p[i]
        end
        ICNp, JCNp, wm2p = ones(Int,sum(cell2 .& phaseMask)), JCN[cell2 .& phaseMask], wm2[phaseMask[cell2]]
        for i in eachindex(ICNp)
            w2_term[ICNp[i],JCNp[i]] += wm2p[i]
        end

        phasesS[Ny,:,n] = ((wc1.*w1_term)./w1 +
                        (wc2.*w2_term)./w2) ./ (wc1+wc2)

        # temp = (wc1*accumarray([ones(sum(cell1 & phaseMask),1), JCN(cell1 & phaseMask)'], wm1(phaseMask(cell1)), [1, Nx])./w1 + ...
        #     wc2*accumarray([ones(sum(cell2 & phaseMask),1), JCN(cell2 & phaseMask)'], wm2(phaseMask(cell2)), [1, Nx])./w2)/...
        #     (wc1+wc2);
        # phasesS(Ny,:,n) = temp;
    end

    ### left edge

    leftEdge = (jcn.==1) .& (icn.>1) .& (icn.<Ny)
    shifted  = (jcn.==1) .& (icn.<Ny-1)

    # cell 1

    cell1 = shifted .& (quad.==4)

    ddx = dx[1]
    ddy = dy[ICN[cell1].-1]
    dxm = xm[cell1] .- x[1]
    dym = ym[cell1] .- y[ICN[cell1]]
    wm1 = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx.*ddy)
    #w1 = accumarray([ICN(cell1)', ones(sum(cell1),1)], wm1, [Ny, 1]);
    w1 = zeros(Float64, (Ny,1))
    ICNp, JCNp = ICN[cell1], ones(Int,sum(cell1))
    for i in eachindex(ICNp)
        w1[ICNp[i],JCNp[i]] += wm1[i]
    end


    # cell 2

    cell2 = leftEdge .& (quad.==1)

    ddx = dx[1]
    ddy = dy[ICN[cell2]]
    dxm = xm[cell2] .- x[1]
    dym = ym[cell2] .- y[ICN[cell2]]
    wm2 = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx.*ddy)
    #w2 = accumarray([ICN(cell2)', ones(sum(cell2),1)], wm2, [Ny, 1]);
    w2 = zeros(Float64, (Ny,1))
    ICNp, JCNp = ICN[cell2], ones(Int,sum(cell2))
    for i in eachindex(ICNp)
        w2[ICNp[i],JCNp[i]] += wm2[i]
    end

    #loop over material properties to interpolate

    Threads.@threads for n = 1:PARAMS.Nphase

        w1_term = zeros(Float64, (Ny,1))
        w2_term = copy(w1_term)

        phaseMask = (phases .== n)

        ICNp, JCNp, wm1p = ICN[cell1 .& phaseMask], ones(Int,sum(cell1 .& phaseMask)), wm1[phaseMask[cell1]]
        for i in eachindex(ICNp)
            w1_term[ICNp[i],JCNp[i]] += wm1p[i]
        end
        ICNp, JCNp, wm2p = ICN[cell2 .& phaseMask], ones(Int,sum(cell2 .& phaseMask)), wm2[phaseMask[cell2]]
        for i in eachindex(ICNp)
            w2_term[ICNp[i],JCNp[i]] += wm2p[i]
        end

        phasesS[:,1,n] = ((wc1.*w1_term)./w1 +
                        (wc2.*w2_term)./w2) ./ (wc1+wc2)

        # temp = (wc1*accumarray([ICN(cell1 & phaseMask)', ones(sum(cell1 & phaseMask),1)], wm1(phaseMask(cell1)), [Ny, 1])./w1 + ...
        #     wc2*accumarray([ICN(cell2 & phaseMask)', ones(sum(cell2 & phaseMask),1)], wm2(phaseMask(cell2)), [Ny, 1])./w2)/...
        #     (wc1+wc2);
        # phasesS(:, 1, n) = temp;
    end

    ### right edge

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
    w1 = zeros(Float64, (Ny,1))
    ICNp, JCNp = ICN[cell1], ones(Int,sum(cell1))
    for i in eachindex(ICNp)
        w1[ICNp[i],JCNp[i]] += wm1[i]
    end

    # cell 2

    cell2 = rightEdge .& (quad.==2)

    ddx = dx[Nx-1]
    ddy = dy[ICN[cell2]]
    dxm = xm[cell2] .- x[Nx-1]
    dym = ym[cell2] .- y[ICN[cell2]]
    wm2 = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx.*ddy)
    #w2 = accumarray([ICN(cell2)', ones(sum(cell2),1)], wm2, [Ny, 1]);
    w2 = zeros(Float64, (Ny,1))
    ICNp, JCNp = ICN[cell2], ones(Int,sum(cell2))
    for i in eachindex(ICNp)
        w2[ICNp[i],JCNp[i]] += wm2[i]
    end

    #loop over material properties to interpolate

    Threads.@threads for n = 1:PARAMS.Nphase

        w1_term = zeros(Float64, (Ny,1))
        w2_term = copy(w1_term)

        phaseMask = (phases .== n)

        ICNp, JCNp, wm1p = ICN[cell1 .& phaseMask], ones(Int,sum(cell1 .& phaseMask)), wm1[phaseMask[cell1]]
        for i in eachindex(ICNp)
            w1_term[ICNp[i],JCNp[i]] += wm1p[i]
        end
        ICNp, JCNp, wm2p = ICN[cell2 .& phaseMask], ones(Int,sum(cell2 .& phaseMask)), wm2[phaseMask[cell2]]
        for i in eachindex(ICNp)
            w2_term[ICNp[i],JCNp[i]] += wm2p[i]
        end

        phasesS[:,Nx,n] = ((wc1.*w1_term)./w1 +
                        (wc2.*w2_term)./w2) ./ (wc1+wc2)

        # temp = (wc1*accumarray([ICN(cell1 & phaseMask)', ones(sum(cell1 & phaseMask),1)], wm1(phaseMask(cell1)), [Ny, 1])./w1 + ...
        #     wc2*accumarray([ICN(cell2 & phaseMask)', ones(sum(cell2 & phaseMask),1)], wm2(phaseMask(cell2)), [Ny, 1])./w2)/...
        #     (wc1+wc2);
        # phasesS(:,Nx, n) = temp;
    end

    ## CORNERS

    # upper left

    upperLeft = (jcn.==1) .& (icn.==1) .& (quad.==1)

    ddx = dx[1]
    ddy = dy[1]
    dxm = xm[upperLeft] .- x[1]
    dym = ym[upperLeft] .- y[1]
    wm  = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx*ddy)
    wco = sum(wm)

    for n = 1:PARAMS.Nphase

        phaseMask = (phases .== n)

        phasesS[1,1,n] = sum(wm[phaseMask[upperLeft]])./wco
    end

    # upper right

    upperRight = (icn.==1) .& (jcn.==(Nx-1)) .& (quad.==2)

    ddx = dx[Nx-1]
    ddy = dy[1]
    dxm = xm[upperRight] .- x[Nx-1]
    dym = ym[upperRight] .- y[1]
    wm  = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx*ddy)
    wco = sum(wm)

    for n = 1:PARAMS.Nphase

        phaseMask = (phases .== n)

        phasesS[1,Nx,n] = sum(wm[phaseMask[upperRight]])./wco
    end

    # lower Right

    lowerRight = (icn.==Ny-1) .& (jcn.==(Nx-1)) .& (quad.==3)

    ddx = dx[Nx-1]
    ddy = dy[Ny-1]
    dxm = xm[lowerRight] .- x[Nx-1]
    dym = ym[lowerRight] .- y[Ny-1]
    wm  = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx*ddy)
    wco = sum(wm)

    for n = 1:PARAMS.Nphase

        phaseMask = (phases .== n)

        phasesS[Ny,Nx,n] = sum(wm[phaseMask[lowerRight]])./wco
    end

    # lower left

    lowerLeft = (icn.==(Ny-1)) .& (jcn.==1) .& (quad.==4)

    ddx = dx[1]
    ddy = dy[Ny-1]
    dxm = xm[lowerLeft] .- x[1]
    dym = ym[lowerLeft] .- y[Ny-1]
    wm  = 1 .- (dxm.*dym .+ (ddx.-dxm).*dym .+ (ddy.-dym).*dxm)./(ddx*ddy)
    wco = sum(wm)

    for n = 1:PARAMS.Nphase

        phaseMask = (phases .== n)

        phasesS[Ny,1,n] = sum(wm[phaseMask[lowerLeft]])./wco
    end

    return phasesS

end
