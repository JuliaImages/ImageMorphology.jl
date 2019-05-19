__precompile__()

module ImageMorphology

using ImageCore

include("convexhull.jl")

include("dilation_and_erosion.jl")
include("thinning.jl")

export
    dilate,
    erode,
    opening,
    closing,
    tophat,
    bothat,
    morphogradient,
    morpholaplace,

    # convexhull.jl
    convexhull,

    # thinning.jl
    thinning,
    GuoAlgo


end # module
