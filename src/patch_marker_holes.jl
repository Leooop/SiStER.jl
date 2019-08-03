# function patch_marker_holes(icn,jcn,quad,Nx,Ny,Mquad,Mquad_crit,xm,ym,x,y,dx,dy,im,ep,idm,Tm,sxxm,sxym,epNH,epsIIm)
#
# seeds new markers in all quadrants where marker density has fallen below
# the threshold
# AND THEN assigns those new markers a marker index, a marker phase and a plastic
# strain based on the properties of the markers surrounding them
# (these particular properties are never passed on the grid, e.g. plastic strain...)
#
# J.-A. Olive, and B.Z. Klein, 2012-2014

function patch_marker_holes(icn,jcn,quad,Nx,Ny,Mquad,Mquad_crit,xm,ym,x,y,dx,dy,phm,ep,idm,Tm,sxxm,sxym,epNH,epsIIm)


    M = length(xm)
    md_crit = Mquad_crit
    mark_per_quad = Mquad
    mdx = minimum(dx)/2
    mdy = minimum(dy)/2

    ############### LOOK FOR EMPTY (no more markers) QUADRANTS

    mp = accumarray(icn, jcn, quad, 1, (Ny-1, Nx-1, 4))
    empty = findall(mp.==0)

    display_message_if_empty = 0
    if display_message_if_empty==1
        if (~isempty(empty))
            #iEmpty,jEmpty,kEmpty = ind2sub(size(mp), empty);
            for n = 1:min(length(empty), 100)
                println("WARNING ! Empty quadrant number ", string(empty[n][3]),  " in cell i = ", string(empty[n][1]), ", j = ", num2str(empty[n][1]))
            end
        end
    end


    ##### LOCATE QUADRANTS WHERE MARKER DENSITY IS BELOW THRESHOLD
    mpCrInd = findall(mp.<md_crit)
    # [iicr, jjcr, qqcr] = ind2sub(size(mp), mpCrInd);
    # iicr = iicr';
    # jjcr = jjcr';
    # qqcr = qqcr';


    # NEED TO SEED NEW MARKERS IN THOSE QUADRANTS
    # SO THAT WE ARE BACK TO THE INITIAL MARKER DENSITY IN THOSE QUADRANTS

    xrsd=Float64[]
    yrsd=Float64[]
    phm_fix=Int[] # marker phase
    ep_fix=Float64[] # plastic strain
    epNH_fix=Float64[] # non-healed plastic strain
    te_fix=Float64[] # temperature
    sxx_fix=Float64[] # stress
    sxy_fix=Float64[] # stress
    sr_fix=Float64[] # strain rate
    #D_fix...

    if ~isempty(mpCrInd) # if there are critical quadrants

        for c in eachindex(mpCrInd) # go through all critical quadrants

            # the upper-left node of the corresponding cell is
            icell = mpCrInd[c][1]
            jcell = mpCrInd[c][2]
            qcell = mpCrInd[c][3]
            # the quadrant area in that cell is
            qsize = 0.25.*dx[jcell].*dy[icell]
            # if the smallest quadrant in the domain (area mdx*mdy) has
            # mark_per_quad markers, then this critical quadrant needs
            Nfix = ceil(qsize/(mdx*mdy))*mark_per_quad
            Nfix = Int(max(Nfix-mp[mpCrInd[c]],1))

            # find the coordinates of the upper-left corner of the quadrant
            if (qcell==1) || (qcell==4)
                xcorn = x[jcell]
            else
                xcorn=x[jcell]+dx[jcell]/2
            end
            if (qcell==1) || (qcell==2)
                ycorn=y[icell]
            else
                ycorn=y[icell]+dy[icell]/2
            end

            # draw random marker location
            xmrr, ymrr = seed_markers_uniformly(xcorn,ycorn,dx[jcell]/2,dy[icell]/2,Nfix)

            #xrsd = [xrsd ; xmrr]
            #yrsd=[yrsd ; ymrr]

            append!(xrsd, xmrr)
            append!(yrsd, ymrr)
    #### NOW THAT THE NEW MARKERS ARE SEEDED,
    ##### ASSIGN PARAMETERS THAT ARE NEVER STORED ON THE EULERIAN GRID

    # the value is assigned based on the average, or max. value of the
    # markers that remain in the corresponding CELL (since there's no grid value to interpolate from)

    # CAREFUL THIS CANNOT WORK IF WE END UP WITH AN EMPTY CELL
    # if that was to happen, let's just draw "im" randomly


        if isempty(ep[(icn.==icell) .& (jcn.==jcell)]) == true

            println("WARNING ! - EMPTY CELL - SOMETHING IS VERY WRONG...")
            phm_fix = 1 .+ floor.(rand(Nfix).*maximum(phm)) # random phase number
            ep_fix = zeros(Float64, Nfix)
            epNH_fix=zeros(Float64, Nfix)
            te_fix=zeros(Float64, Nfix)
            sxx_fix=zeros(Float64, Nfix)
            sxy_fix=zeros(Float64, Nfix)
            sr_fix=zeros(Float64, Nfix)
            #D_fix= ...

        else


            # assign the average phase of the markers that are left in the cell
            phase_fix = round(mode(phm[(icn.==icell) .& (jcn.==jcell)]))
            # assign the greatest plastic strain of the markers that are left in the cell
            strain_fix = maximum(ep[icn.==icell .& jcn.==jcell])
            strainNH_fix = maximum(epNH[icn.==icell .& jcn.==jcell])
            # assign the average temperature of the markers that are left in
            # the cell
            temp_fix = mean(Tm[icn.==icell .& jcn.==jcell])
            # reassign mean stress / strain rate of markers left in cell
            stress_xx_fix = mean(sxxm[icn.==icell .& jcn.==jcell])
            stress_xy_fix = mean(sxym[icn.==icell .& jcn.==jcell])
            strainrate_fix = mean(epsIIm[icn.==icell .& jcn.==jcell])
            #damage_fix = ... similar to Ep (max)

        end

        append!(phm_fix, phase_fix.*ones(Int,Nfix))
        append!(ep_fix, strain_fix.*ones(Float64,Nfix))
        append!(epNH_fix, strainNH_fix.*ones(Float64,Nfix))
        append!(te_fix, temp_fix.*ones(Float64,Nfix))
        append!(sxx_fix, stress_xx_fix.*ones(Float64,Nfix))
        append!(sxy_fix, stress_xy_fix.*ones(Float64,Nfix))
        append!(sr_fix, strainrate_fix.*ones(Float64,Nfix))

        #D_fix = ...

        end



    # NOW ASSIGN PROPERTIES TO THOSE MARKERS
    #Npatch = length(xrsd)
    #index_fix = (maximum(idm) + 1):(maximum(idm) + Npatch)

    #Ifix = (M+1):(M+Npatch) # total number of markers added to fix holes in critical quadrants
    append!(xm,xrsd)
    append!(ym,yrsd)
    append!(im,im_fix)
    append!(ep,ep_fix)
    append!(epNH,epNH_fix)
    append!(idm,index_fix)
    append!(Tm,te_fix)
    append!(sxxm,sxx_fix)
    append!(sxym,sxy_fix)
    append!(epsIIm,sr_fix)
    #Dm = ...

    # uncomment to display number of added markers
    #fprintf('\n#d#s#d#s\n', length(Ifix), ' markers added in ', length(iicr), ' cell quadrants.')

    else

    Ifix = 0

    end

    return xm, ym, im, Ifix, mp, ep, idm, Tm, sxxm, sxym, epNH, epsIIm

end
