module SiStER

#--------------
# EXPORTATION #
#--------------

#----------------
# PACKAGES USED #
#----------------
# standard :
using
    Random: randperm,
    Statistics: mean

# third-party :
using
    LazyJSON: PropertyDict,
    Interpolations

#--------
# FILES #
#--------

include("main.jl")
end # module
