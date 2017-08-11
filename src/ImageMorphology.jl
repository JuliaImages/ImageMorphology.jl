__precompile__(true)

module ImageMorphology

using ImageCore
using ImageMetadata

include("morphological_operations.jl")

export
    dilate,
    erode,
    opening,
    closing,
    tophat,
    bothat,
    morphogradient,
    morpholaplace

end # module
