__precompile__()

module ImageMorphology

using ImageCore
using Images

abstract type SkeletonizationAlgo end
struct MedialAxis <: SkeletonizationAlgo end


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
    MedialAxis,
    skeletonize,
    thinning,
    GuoAlgo


end # module
