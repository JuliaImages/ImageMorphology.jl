module ImageMorphology

using DataStructures: Queue, enqueue!, dequeue!
using ImageCore
using ImageCore: GenericGrayImage
using OffsetArrays
using OffsetArrays: centered
using LinearAlgebra
using TiledIteration: EdgeIterator, SplitAxis, SplitAxes
using Requires
using LoopVectorization

const _docstring_se = """
`se` is the structuring element that defines the neighborhood of the image. See
[`strel`](@ref) for more details. If `se` is not specified, then it will use the
[`strel_box`](@ref) with an extra keyword `dims` to control the dimensions to filter,
and half-size `r` to control the diamond size.
"""
include("StructuringElements/StructuringElements.jl")
using .StructuringElements

include("convexhull.jl")
include("connected.jl")
include("clearborder.jl")
include("extreme_filter.jl")
include("ops/dilate.jl")
include("ops/erode.jl")
include("ops/closing.jl")
include("ops/opening.jl")
include("ops/tophat.jl")
include("ops/bothat.jl")
include("ops/mgradient.jl")
include("ops/mlaplacian.jl")
include("ops/mreconstruct.jl")
include("ops/underbuild.jl")
include("ops/overbuild.jl")
include("isboundary.jl")
include("thinning.jl")
include("imfill.jl")
include("maxtree.jl")
include("leveling.jl")

include("feature_transform.jl")
include("utils.jl")
using .FeatureTransform

include("deprecations.jl")

export
    # structuring_element.jl
    centered, # from OffsetArrays
    strel,
    strel_chain,
    strel_product,
    strel_type,
    strel_size,
    strel_diamond,
    strel_box,

    # operations
    dilate,
    dilate!,
    erode,
    erode!,
    extreme_filter,
    extreme_filter!,
    opening,
    opening!,
    closing,
    closing!,
    tophat,
    tophat!,
    bothat,
    bothat!,
    mgradient,
    mgradient!,
    mlaplacian,
    mlaplacian!,
    mreconstruct,
    mreconstruct!,
    underbuild,
    underbuild!,
    overbuild,
    overbuild!,

    # connected.jl
    label_components,
    label_components!,
    component_boxes,
    component_lengths,
    component_indices,
    component_subscripts, # deprecated (v0.4)
    component_centroids,

    # convexhull.jl
    convexhull,

    # isboundary.jl
    isboundary,
    isboundary!,

    # thinning.jl
    thinning,
    GuoAlgo,
    imfill,

    # maxtree.jl
    MaxTree,
    areas,
    boundingboxes,
    diameters,
    area_opening,
    area_opening!,
    area_closing,
    area_closing!,
    diameter_opening,
    diameter_opening!,
    diameter_closing,
    diameter_closing!,
    local_maxima!,
    local_maxima,
    local_minima!,
    local_minima,

    #feature_transform.jl
    feature_transform,
    distance_transform,

    #leveling
    low_leveling,
    low_leveling!,
    high_leveling,
    high_leveling!,
    leveling,
    leveling!, clearborder

function __init__()
    @require ImageMetadata = "bc367c6b-8a6b-528e-b4bd-a4b897500b49" begin
        # morphological operations for ImageMeta
        function dilate(img::ImageMetadata.ImageMeta; kwargs...)
            out = dilate!(similar(ImageMetadata.arraydata(img)), img; kwargs...)
            return ImageMetadata.shareproperties(img, out)
        end
        function erode(img::ImageMetadata.ImageMeta; kwargs...)
            out = erode!(similar(ImageMetadata.arraydata(img)), img; kwargs...)
            return ImageMetadata.shareproperties(img, out)
        end
    end
end

end # module
