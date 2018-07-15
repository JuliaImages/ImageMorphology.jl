__precompile__(true)

module ImageMorphology

using ImageCore

include("morphological_operations.jl")
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

    thinning,
    GuoAlgo
    
end # module
