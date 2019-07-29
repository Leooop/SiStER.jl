using Random: randperm
using Statistics: mean
using SparseArrays

import Base.Threads.@spawn
import LinearAlgebra.norm

# third-party :
using LazyJSON: PropertyDict
using Interpolations

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
