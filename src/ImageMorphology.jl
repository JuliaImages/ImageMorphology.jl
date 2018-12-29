__precompile__()

module ImageMorphology

using ImageCore

include("dilation_and_erosion.jl")
include("thinning.jl")
include("imfill.jl")

export
    dilate,
    erode,
    opening,
    closing,
    tophat,
    bothat,
    morphogradient,
    morpholaplace,

    thinning,
    GuoAlgo
    imfill


end # module
