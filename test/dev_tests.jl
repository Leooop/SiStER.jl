function test()

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
end
@btime test()


function loop1()
    phasesS = Array{Float64,3}(undef,Ny,Nx,PARAMS.Nphase)
    for n = 1:PARAMS.Nphase
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
    return phasesS
end

@everywhere function inner_func(n,Nx,Ny,ICN,JCN,wm1,wm2,wm3,wm4,w1,w2,w3,w4)
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

    return ((wc1.*w1_term)./w1 +
                    (wc2.*w2_term)./w2 +
                    (wc3.*w3_term)./w3 +
                    (wc4.*w4_term)./w4 )./ (wc1+wc2+wc4+wc4)
end

function loopn()
    begin
        phases_vec = collect(1:PARAMS.Nphase)
        interm = pmap(n->inner_func(n), 1:PARAMS.Nphase, distributed = false)
        cat(interm[1],interm[2],interm[3],dims = 3)
    end
end

function loopn2()
    interm = Vector{Array{Float64,2}}(undef,PARAMS.Nphase)
    @sync begin
        @async for n = 1:PARAMS.Nphase
                    f = @spawn inner_func(n,Nx,Ny,ICN,JCN,wm1,wm2,wm3,wm4,w1,w2,w3,w4)
                    interm[n] = fetch(f)
                end
    end
    return cat(interm[1],interm[2],interm[3],dims = 3)
end

using SharedArrays
function loop3n()
    phasesS = SharedArray{Float64,3}((Ny,Nx,PARAMS.Nphase))
    @distributed for n = 1:PARAMS.Nphase
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
    return phasesS
end

p1 = @benchmark loop1()
pn = @benchmark loopn()
p2n = @benchmark loop2n()
p3n = @benchmark loop3n()


using SharedArrays
using Distributed

n = 100000

function testn(n)
    A = Vector{Float64}(undef,n)
    @time Threads.@threads for i = 1:n
        A[i] = sqrt(i/pi)
    end
end

function test1(n)
    A = Vector{Float64}(undef,n)
    @time for i = 1:n
        A[i] = sqrt(i/pi)
    end
end

testn(n)
test1(n)
