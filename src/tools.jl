####Tools functions

function sub2ind(lin,i,j)
    lin[i,j]
end

using SharedArrays
using Distributed

n = 10000000
As = SharedArray{Float64}(n)
@time @distributed for i = 1:n
    A[i] = i
end

A = Vector{Float64}(undef,n)
@time for i = 1:n
    A[i] = i
end
