# main function
#
# Simple Stokes solver with Exotic Rheologies
#
# Main routine doing initialization, time loop and outputs
#
#
# J.-A. Olive, B.Z. Klein, E. Mittelstaedt, M. Behn, G. Ito, S. Howell
# jaolive <at> ldeo.columbia.edu
# March 2011 - April 2017


function main(InpFil::String)

    # INITIALIZATION

    # Input File: loads parameter values, model geometry, boundary conditions
    include(InpFil)

    # construct grid and initialize marker / node arrays
    include("initialize.jl")

    # BEGIN TIME LOOP #########################################################
    time_sim = 0

    for t=1:Nt # time loop

        println("STARTING ITERATION: ", string(t), " out of ", string(Nt))

        # update time
        time_sim += dt_m

        # Here we prepare nodal arrays to feed the Stokes solver
        include("material_props_on_nodes.jl")

        ### SOLVE STOKES WITH NON-LINEAR RHEOLOGY HERE
        include("flow_solve.jl")

        # GET STRAIN RATE FROM CURRENT SOLUTION
        epsIIm = interp_shear_nodes_to_markers(epsII_s,x,y,xm,ym,icn,jcn)

        # USE STRAIN RATE TO UPDATE STRESSES ON MARKERS
        include("update_marker_stresses.jl")

        # BUILD UP PLASTIC STRAIN IN YIELDING AREAS IF PLASTICITY IS ACTIVATED
        if (PARAMS.YNPlas==1)
            update_ep
        end

        # OUTPUT VARIABLES OF INTEREST (prior to rotation & advection)
        if (mod(t,dt_out)==0 && dt_out>0) || t==1 || t==Nt # SAVING SELECTED OUTPUT
            println("SAVING SELECTED VARIABLES TO OUTPUT FILE")
            filename=string(t)
            etam = interp_shear_nodes_to_markers(etas,x,y,xm,ym,icn,jcn) # to visualize viscosity on markers
            @save "$filename.jld2" X Y vx vy p time xm ym etam rhom BC etan Tm im idm epsIIm sxxm sxym ep epNH icn jcn qd
        end

        # SET ADVECTION TIME STEP BASED ON CURRENT FLOW SOLUTION
        dt_m = set_timestep(dx,dy,vx,vy,PARAMS)

        # ROTATE ELASTIC STRESSES IN CURRENT FLOW FIELD
        if (PARAMS.YNElast==1)
            rotate_stresses
        end

        # EVOLVE TEMPERATURE FIELD THROUGH DIFFUSION
        if PARAMS.Tsolve==1
            SiStER_thermal_update;
        end

        # MARKER ADVECTION, REMOVAL, AND ADDITION #############################
        SiStER_move_remove_and_reseed_markers;
        # advect markers in current flow field
        # remove markers if necessary
        # add markers if necessary
        SiStER_update_topography_markers
        # here we do the same for the marker chain that keeps track of topography
        #######################################################################

        println('---------------')
        println(['END OF ITERATION: ' num2str(t) ' out of ' num2str(Nt) ' - SIMULATION TIME: ' num2str(time/365.25/24/3600/1000) ' kyrs.'])
        println('--------------------------------')
        println('--------------------------------')


    end

end
println('FIN')
