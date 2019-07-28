####Tools functions

sub2ind(lin,i,j) = lin[i,j]

sub2ind(lin,cart) = lin[cart]

function accumarray(inds::Vector, val::Vector, size::Int)
    acc = zeros(Float64, size)
    for i in eachindex(inds)
        acc[inds[i]] += val[i]
    end
    return acc
end

function accumarray(inds1::Vector, inds2::Vector, val::Vector, size::Tuple)
    acc = zeros(Float64, size)
    for i in eachindex(inds1)
        acc[inds1[i],inds2[i]] += val[i]
    end
    return acc
end
