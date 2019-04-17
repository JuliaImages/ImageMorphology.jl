__precompile__()

module ImageMorphology

using ImageCore
using Images

abstract type AbstractSkeletonizationAlgorithm end
struct MedialAxisTransform <: AbstractSkeletonizationAlgorithm end


include("dilation_and_erosion.jl")
include("thinning.jl")
include("skeletonization.jl")
export
    dilate,
    erode,
    opening,
    closing,
    tophat,
    bothat,
    morphogradient,
    morpholaplace,
    MedialAxisTransform,
    skeletonize,
    thinning,
    GuoAlgo


end # module
