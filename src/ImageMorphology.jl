__precompile__()

module ImageMorphology

using ImageCore
using Base.Cartesian # TODO: delete this

include("convexhull.jl")
include("connected.jl")

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

    # connected.jl
    label_components,
    label_components!,

    # convexhull.jl
    convexhull,

    # thinning.jl
    thinning,
    GuoAlgo


end # module
