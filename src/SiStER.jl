module SiStER

#--------------
# EXPORTATION #
#--------------

#----------------
# PACKAGES USED #
#----------------
# standard :
using Random: randperm
using Statistics: mean

# third-party :
using LazyJSON: PropertyDict
using Interpolations

#--------
# FILES #
#--------

include("functions_load.jl")
include("main.jl")
end # module
