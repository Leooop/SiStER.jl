# [km, cpm]=SiStER_get_thermal_properties(im,MAT)
# obtain thermal conductivity and heat capacity

function get_thermal_properties(phm,MAT)

    km = zeros(Float64,size(phm))
    cpm = copy(km)

    types = unique(phm)
    for i in eachindex(types)
        imInd = (phm .== types[i])
        km[imInd] .= MAT[types[i]].k
        cpm[imInd] .= MAT[types[i]].cp
    end

    return km, cpm

end
