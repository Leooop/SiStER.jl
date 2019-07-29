module SiStER

#--------------
# EXPORTATION #
#--------------

##----------------
# PACKAGES USED #
#----------------
# standard :
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
include("main.jl")
end # module
