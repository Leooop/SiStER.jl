# [quad,icn,jcn] = SiStER_locate_markers_in_grid(xm,ym,x,y,dx,dy)
# Tells a marker which cell (and which quadrant of that cell) it belongs to.
    # icn,jcn are the indexes of the upper-left shear node of the cell
    # that a given marker is currently in
    # quad is the quadrant of the cell that contains the marker
    # (quad = 1 means bottom-right, then numbered clockwise)
    # sped up by B. Klein in Fall 2016 by using interp1 function

function locate_markers_in_grid(xm::Vector{Float64},ym::Vector{Float64},x::Vector{Float64},y::Vector{Float64},dx,dy)

    ## Determine Location of Markers and Quadrant of Element
    M=length(xm)
    icn = zeros(Int64,M)
    jcn = copy(icn)
    quad = copy(icn) # quadrant 1 = bottom-right, numbered clockwise

    # Create interpolation object for the index with nearnest (Constant()) parametrization
    itpX = interpolate((x,), 1:length(x), Gridded(Constant()))
    itpY = interpolate((y,), 1:length(y), Gridded(Constant()))

    # Allow extrapolation with values outside the boundaries being equal to the boundaries values
    itpX = extrapolate(itpX, Flat())
    itpY = extrapolate(itpY, Flat())

    # Interpolate to {xm,ym}
    Ix = Int.(itpX(xm))
    Iy = Int.(itpY(ym))

    # Old strategy :
    # [~, Ix] = min(abs(bsxfun(@minus, xm, x')));
    # [~, Iy] = min(abs(bsxfun(@minus, ym, y')));

    jcn[xm.>x[Ix]]  = Ix[xm.>x[Ix]]
    jcn[xm.<=x[Ix]] = Ix[xm.<=x[Ix]].-1

    icn[ym.>y[Iy]]  = Iy[ym.>y[Iy]]
    icn[ym.<=y[Iy]] = Iy[ym.<=y[Iy]].-1

    jcn[jcn.==0] .= 1
    jcn[jcn.>length(dx)] .= length(dx)

    icn[icn.==0] .= 1
    icn[icn.>length(dy)] .= length(dy)

    disx = abs.((xm.-x[jcn])./dx[jcn])
    disy = abs.((ym.-y[icn])./dy[icn])


    xRIGHT = disx .> 0.5
    yUP = disy .> 0.5

    quad[xRIGHT .& yUP] .= 3
    quad[xRIGHT .& .~yUP]  .= 2
    quad[.~xRIGHT .& yUP]  .= 4
    quad[.~xRIGHT .& .~yUP]   .= 1

    return quad, icn, jcn

end

# Visualization :
# scatter(X,Y)
# inds = (icn .>= 0) .& (jcn .>= 0) .& (qd .== 1)
# scatter!(xm[inds],ym[inds],ms = 0.5)
