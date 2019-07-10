function get_density(phm,Tm,MAT)
    # [rho]=SiStER_get_density(im,Tm,MAT)
    # obtain density from temperature and material identity
    #
    # edited by B. Klein, 9/20/2016 to remove struct indexing concatenating

    T0=0.0
    rhom = zeros(Float64,size(phm))

    types = unique(phm)
    for i in eachindex(types)
        logical = (phm .== types[i])
        rho0 = MAT[types[i]].rho0
        alpha = MAT[types[i]].alpha
        rhom[logical] = rho0.*(1 .-alpha.*(Tm[logical] .-T0))
    end

    return rhom
    # OLD WAY (slow)
    #T0=0;
    #rhom = [MAT(im).rho0].*(1-[MAT(im).alpha].*(Tm-T0));
end
