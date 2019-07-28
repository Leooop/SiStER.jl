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

import Base.Threads.@spawn

# third-party :
using LazyJSON: PropertyDict
using Interpolations

##--------
# FILES #
#--------

include("functions_load.jl")
include("main.jl")
end # module
