# [ys]=SiStER_get_elastic_moduli(im,MAT)
# obtain shear modulus from material properties
# B. Klein 9/16

function get_elastic_moduli(phm,MAT)

    Gm = zeros(Float64,size(phm))

    types = unique(phm)
    for i = 1:length(types)
        imInd = (phm .== types[i])
        Gm[imInd] .= MAT[types[i]].G
    end
    return Gm

end

# faster than:
# Gm=[MAT(im).G];
