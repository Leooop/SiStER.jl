####Tools functions

function sub2ind(lin,i,j)
    lin[i,j]
end

# using SharedArrays
# using Distributed
#
# n = 100000000
#
# function testn(n)
#     A = Vector{Float64}(undef,n)
#     @time Threads.@threads for i = 1:n
#         A[i] = sqrt(i/pi)
#     end
# end
#
# function test1(n)
#     A = Vector{Float64}(undef,n)
#     @time for i = 1:n
#         A[i] = sqrt(i/pi)
#     end
# end
#
# testn(n)
# test1(n)
