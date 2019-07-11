# [fric]=SiStER_get_friction(im,ep,MAT)
# compute friction on markers based on plastic strain, like cohesion
# J.-A. Olive 4/2017

function get_friction(phm,ep,MAT)

    fric = Vector{Float64}(undef,size(phm))

    types = unique(phm)
    for i in eachindex(types)
        imInd = (phm .== types[i])
        mumax=MAT[types[i]].mu

        if in(:mumin,keys(MAT[1]))
            mumin=MAT[types[i]].mumin
        else
            mumin=MAT[types[i]].mu
        end

        epscrit=MAT[types[i]].ecrit
        muminvec = fill(mumin,count(imInd))
        # get cohesion
        fric[imInd]=maximum((mumax.+(mumin.-mumax).*ep[imInd]./epscrit,muminvec))
    end

    return fric

end

## OLD (slow due to brackets)
# mumax=[MAT(im).mu];
#
#
# if isfield(MAT,'mumin')==1
#     mumin=[MAT(im).mumin];
# else
#     mumin=[MAT(im).mu];
# end
#
# epscrit=[MAT(im).ecrit];
#
# # get cohesion
# fric=max(mumax+(mumin-mumax).*ep./epscrit,mumin);
#
# return
