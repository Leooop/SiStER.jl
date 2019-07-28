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

function acump()
    n2interp = Vector{Matrix{Float64}}(undef,1)
    for vn = 1:1

        w1_term = Threads.@spawn accumarray(ICN[cell1], JCN[cell1],
            args[vn][cell1].*wm1, (Ny,Nx))
        w2_term = Threads.@spawn accumarray(ICN[cell2], JCN[cell2],
            args[vn][cell2].*wm2, (Ny,Nx))
        w3_term = Threads.@spawn accumarray(ICN[cell3], JCN[cell3],
            args[vn][cell3].*wm3, (Ny,Nx))
        w4_term = Threads.@spawn accumarray(ICN[cell4], JCN[cell4],
            args[vn][cell4].*wm4, (Ny,Nx))

        # wait(t1)
        # wait(t2)
        # wait(t3)
        # wait(t4)

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
    return n2interp
end

function acum()
    n2interp = Vector{Matrix{Float64}}(undef,1)
    for vn = 1:1

        w1_term = accumarray(ICN[cell1], JCN[cell1], args[vn][cell1].*wm1, (Ny,Nx))
        w2_term = accumarray(ICN[cell2], JCN[cell2], args[vn][cell2].*wm2, (Ny,Nx))
        w3_term = accumarray(ICN[cell3], JCN[cell3], args[vn][cell3].*wm3, (Ny,Nx))
        w4_term = accumarray(ICN[cell4], JCN[cell4], args[vn][cell4].*wm4, (Ny,Nx))

        n2interp[vn] = ((wc1.*w1_term)./w1 +
                            (wc2.*w2_term)./w2 +
                            (wc3.*w3_term)./w3 +
                            (wc4.*w4_term)./w4 )./ (wc1+wc2+wc4+wc4)

        # n2interp(vn).data = (wc1*accumarray([ICN(cell1)', JCN(cell1)'], args{vn}(cell1).*wm1)./w1 + ...
        #     wc2*accumarray([ICN(cell2)', JCN(cell2)'], args{vn}(cell2).*wm2)./w2 + ...
        #     wc3*accumarray([ICN(cell3)', JCN(cell3)'], args{vn}(cell3).*wm3)./w3 + ...
        #     wc4*accumarray([ICN(cell4)', JCN(cell4)'], args{vn}(cell4).*wm4)./w4)./...
        #     (wc1+wc2+wc4+wc4);
    end
    return n2interp
end


p1 = @benchmark acum()
ps = @benchmark acump()



using SharedArrays
using Distributed

n = 100000

function testn(n)
    A = Vector{Float64}(undef,n)
    Threads.@threads for i = 1:n
        A[i] = sqrt(i/pi)
    end
    return A
end

function testn1(n)
    A = Vector{Float64}(undef,n)
    Threads.@spawn for i = 1:n
        A[i] = sqrt(i/pi)
    end
    return A
end

function test1(n)
    A = Vector{Float64}(undef,n)
    for i = 1:n
        A[i] = sqrt(i/pi)
    end
    return A
end

function testn3(n)
    A = 0
    p = Threads.@spawn for i = 1:n
        A += i
    end
    #wait(p)
    return A
end

function test3(n)
    A = 0
    for i = 1:n
        A += i
    end
    return A
end

@time test1(n)
@time testn(n)
@time testn1(n)
@time testn3(n)
@time test3(n)

test = [1 2 ; 3 4]
fv() = @view test[:,1]
f() = test[:,1]

@benchmark fv()
@benchmark f()

function acum()
    ICNp, JCNp = ICN[cell1], JCN[cell1]
    w1 = zeros(Float64, (Ny,Nx))
    for i in eachindex(ICNp)
        w1[ICNp[i],JCNp[i]] += wm1[i]
    end
end
@benchmark acum()

@benchmark accumarray(ICN[cell1], JCN[cell1], wm1, (Ny,Nx))
