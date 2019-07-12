# get_ductile_rheology_on_nodes_from_strain_rate(MAT,PARAMS,T,epsII,phase_node)
# computes ductile rheology given phase numbers around a node or cell
# (phase_node is either phase_n or phase_s)
#
# G.Ito 8/2016 # JAO 9/2016 fixed parentheses in exp(Ea/(nRT))
# B.Z. Klein 07/17: can now handle nodes surrounded by any number of phases
# that are NOT consecutively numbered (now looping through all phases)
# thanks to new phase interp functions.

function get_ductile_rheology_on_nodes_from_strain_rate(MAT,PARAMS,T,epsII,phase_node)

    pre_diff_vec = [Float64(MAT[i].pre_diff) for i = 1:PARAMS.Nphase]
    ndiff_vec = [Int(MAT[i].ndiff) for i = 1:PARAMS.Nphase]
    Ediff_vec = [Float64(MAT[i].Ediff) for i = 1:PARAMS.Nphase]
    pre_disc_vec = [Float64(MAT[i].pre_disc) for i = 1:PARAMS.Nphase]
    ndisc_vec = [Int(MAT[i].ndisc) for i = 1:PARAMS.Nphase]
    Edisc_vec = [Float64(MAT[i].Edisc) for i = 1:PARAMS.Nphase]

    pre_diff = permutedims(repeat(pre_diff_vec, outer = [1,size(T)...]), [2,3,1])
    ndiff = permutedims(repeat(ndiff_vec, outer = [1,size(T)...]), [2,3,1])
    Ediff = permutedims(repeat(Ediff_vec, outer = [1,size(T)...]), [2,3,1])
    pre_disc = permutedims(repeat(pre_disc_vec, outer = [1,size(T)...]), [2,3,1])
    ndisc = permutedims(repeat(ndisc_vec, outer = [1,size(T)...]), [2,3,1])
    Edisc = permutedims(repeat(Edisc_vec, outer = [1,size(T)...]), [2,3,1])

    epsII = repeat(epsII, outer = [1, 1, PARAMS.Nphase])
    Tr = repeat(T, outer = [1, 1, PARAMS.Nphase])


    # TO BE REPLACED BY SiStER_flow_law_function SOON
    #eta_diff = pre_diff.^(-1./ndiff).*epsII.^((1-ndiff)./ndiff).* ...
    #    exp((Ediff)./(ndiff.*PARAMS.R.*(T+273.15)));
    #eta_disc = pre_disc.^(-1./ndisc).*epsII.^((1-ndisc)./ndisc).* ...
    #    exp((Edisc)./(ndisc.*PARAMS.R.*(T+273.15)));

    eta_diff = flow_law_function("from_strain_rate",pre_diff,Ediff,ndiff,PARAMS.R,Tr,epsII,zeros(Float64, size(epsII)),PARAMS)
    eta_disc = flow_law_function("from_strain_rate",pre_disc,Edisc,ndisc,PARAMS.R,Tr,epsII,zeros(Float64, size(epsII)),PARAMS)


    # linearly average between viscosity of each phase type
    eta_diffAVG = sum(eta_diff.*phase_node, dims = 3)
    eta_discAVG = sum(eta_disc.*phase_node, dims = 3)

    eta = (1 ./ eta_diffAVG .+ 1 ./eta_discAVG).^(-1)
    eta[eta.<PARAMS.etamin] .= PARAMS.etamin

    return eta

end
