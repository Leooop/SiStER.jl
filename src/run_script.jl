using Random: randperm
using Statistics: mean
using SparseArrays

import Base.Threads.@spawn
import LinearAlgebra.norm

# third-party :
using LazyJSON: PropertyDict
using Interpolations
using JLD2, FileIO # To save variables.

##--------
# FILES #
#--------

include("functions_load.jl")

include("../Input_File_continental_rift.jl")
include("initialize.jl")

# BEGIN TIME LOOP #########################################################
time_sim = 0

t = 1

time_sim += dt_m

# Here we prepare nodal arrays to feed the Stokes solver
include("material_props_on_nodes_parallel.jl")
include("flow_solve.jl")


epsIIm = interp_shear_nodes_to_markers(epsII_s,x,y,xm,ym,icn,jcn)

# USE STRAIN RATE TO UPDATE STRESSES ON MARKERS
include("update_marker_stresses.jl")

# BUILD UP PLASTIC STRAIN IN YIELDING AREAS IF PLASTICITY IS ACTIVATED
if (PARAMS.YNPlas==1)
    include("update_ep.jl")
end

# SET ADVECTION TIME STEP BASED ON CURRENT FLOW SOLUTION
dt_m = set_timestep(dx,dy,vx,vy,PARAMS)

# ROTATE ELASTIC STRESSES IN CURRENT FLOW FIELD
if (PARAMS.YNElast==1)
    include("rotate_stresses.jl")
end

if PARAMS.Tsolve==1
    include("thermal_update.jl")
end
