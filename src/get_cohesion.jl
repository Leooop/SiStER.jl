# get_cohesion(im,ep,MAT)
# compute cohesion on markers based on ep
# G.Ito 8/2016
# sped up by B. Klein 9/2016

function get_cohesion(phm,ep,MAT)

    cohes = Vector{Float64}(undef,size(phm))

    types = unique(phm)
    for i in eachindex(types)
        imInd = (phm .== types[i])
        Cmax=MAT[types[i]].Cmax
        Cmin=MAT[types[i]].Cmin
        Cminvec = fill(Cmin,count(imInd))
        epscrit=MAT[types[i]].ecrit

        # get cohesion
        cohes[imInd]=maximum((Cmax.+(Cmin.-Cmax).*ep[imInd]./epscrit,Cminvec))
    end

    return cohes

end

## OLD VERSION the bracket approach is really slow
#
# Cmax=[MAT(im).Cmax];
# Cmin=[MAT(im).Cmin];
# epscrit=[MAT(im).ecrit];
#
# # get cohesion
# cohes=max(Cmax+(Cmin-Cmax).*ep./epscrit,Cmin);
#
# return
