__precompile__(true)

module ImageMorphology

using ImageCore
using ImageFiltering
include("morphological_operations.jl")

export
    dilate,
    erode,
    opening,
    closing,
    tophat,
    bothat,
    morphogradient,
    morpholaplace,
    thinning
    
end # module
