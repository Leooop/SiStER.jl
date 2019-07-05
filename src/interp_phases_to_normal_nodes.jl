# interp_phases_to_normal_nodes(xm,ym,icn,jcn,x,y,phases, PARAMS.Nphase)
# B.Z. Klein, July 2017, an interp function specific to phase (for normal nodes), to enable
# exact mixing of several phases
# this function is slower than on matlab. Should be improved.

function interp_phases_to_normal_nodes(xm,ym,icn,jcn,x,y,phases,PARAMS)

    Nx=length(x)
    Ny=length(y)
    dx=diff(x)
    dy=diff(y)
    Nicn = length(icn)

    lin = LinearIndices((Ny-1,Nx-1))
    INDEX = sub2ind.(Ref(lin),icn,jcn)

    #original : AcCell = bsxfun(@times, dy', dx);
    AcCell = dy*dx'

    xN = x[1:Nx-1] .+ dx./2
    yN = y[1:Ny-1] .+ dy./2

    # original : [XN, YN] = meshgrid(xN, yN);
    XN = repeat(xN',outer = size(yN))
    YN = repeat(yN,outer = size(xN'))

    AMvec = abs.((xm .- XN[INDEX]).*(ym .- YN[INDEX]))
    WMvec = (AcCell[INDEX] .- AMvec)./AcCell[INDEX]

    phaseWeights=zeros(Float64, (Ny,Nx,PARAMS.Nphase))

    for n = 1:PARAMS.Nphase

        phaseMask = phases.==n
        #phaseWeights[2:Ny,2:Nx,n] = accumarray([icn(phaseMask) jcn(phaseMask)], WMvec(phaseMask), size(AcCell), @sum);
        icnp, jcnp, WMvecp = icn[phaseMask], jcn[phaseMask], WMvec[phaseMask]
        for i in eachindex(icnp)
            phaseWeights[icnp[i]+1,jcnp[i]+1,n] += WMvecp[i]
        end
    end

    sumWeights = repeat(sum(phaseWeights, dims = 3), outer = (1,1,PARAMS.Nphase))
    phaseWeights = phaseWeights./sumWeights

    return phaseWeights

end
