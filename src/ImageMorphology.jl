__precompile__()

module ImageMorphology

using ImageCore

include("dilation_and_erosion.jl")

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
