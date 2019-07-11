# [n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,varargin)
# interpolates properties (in the order of input) from markers to normal nodes
#
# First cut - J.A. Olive, March 2011
# Modified by E. Mittelstaedt, April 2011 to allow multiple inputs
# Modified by B.Z. Klein, Spring 2014 for speedup


function interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,varargin...)

    Nx=length(x)
    Ny=length(y)
    dx=diff(x)
    dy=diff(y)

    # check for number of properties to interpolate
    numV = length(varargin)
    n2interp = Vector{Matrix{Float64}}(undef,numV)

    lin = LinearIndices((Ny-1,Nx-1))
    INDEX = sub2ind.(Ref(lin),icn,jcn)

    AcCell = dy*dx'

    xN = x[1:Nx-1] .+ dx./2
    yN = y[1:Ny-1] .+ dy./2

    # original : [XN, YN] = meshgrid(xN, yN);
    XN = repeat(xN',outer = size(yN))
    YN = repeat(yN,outer = size(xN'))

    AMvec = abs.((xm .- XN[INDEX]).*(ym .- YN[INDEX]))
    WMvec = (AcCell[INDEX] .- AMvec)./AcCell[INDEX]

    #w = accumarray([icn' jcn'], WMvec', [], @sum);
    w = zeros(Float64,Ny-1,Nx-1)
    for i in eachindex(icn)
        w[icn[i],jcn[i]] += WMvec[i]
    end

    for vn = 1:numV
        VecData = varargin[vn].*WMvec

        n2interp[vn] = zeros(Float64,Ny,Nx)
        for i in eachindex(icn)
            n2interp[vn][icn[i]+1,jcn[i]+1] += VecData[i]./w[icn[i],jcn[i]]
        end

        #n2interp[vn][2:Ny,2:Nx] = accumarray([icn' jcn'], VecData', [], @sum)./w;
    end

    return [PropertyDict(Dict(:data => n2interp[i])) for i in eachindex(n2interp)]
end
