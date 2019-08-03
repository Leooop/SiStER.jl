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

function accumarray(inds::Tuple, val::Vector)
    acc = zeros(Float64, maximum(inds[1]), maximum(inds[2]))
    for i in eachindex(inds[1])
        acc[inds[1][i],inds[2][i]] += val[i]
    end
    return acc
end

function accumarray(inds1::Vector, inds2::Vector, inds3::Vector, val::Vector, size::Tuple)
    acc = zeros(Float64, size)
    for i in eachindex(inds1)
        acc[inds1[i],inds2[i],inds3[i]] += val[i]
    end
    return acc
end

function accumarray(inds1::Vector, inds2::Vector, inds3::Vector, val::Real, size::Tuple)
    acc = zeros(Float64, size)
    for i in eachindex(inds1)
        acc[inds1[i],inds2[i],inds3[i]] += val
    end
    return acc
end
