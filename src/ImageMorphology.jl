module ImageMorphology

using ImageCore
using ImageCore: GenericGrayImage
using ImageCore.OffsetArrays
using ImageCore.OffsetArrays: centered
using LinearAlgebra
using TiledIteration: EdgeIterator, SplitAxis, SplitAxes
using Requires

include("structuring_element.jl")
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
include("ops/morpholaplace.jl")
include("isboundary.jl")
include("thinning.jl")
include("imfill.jl")
include("maxtree.jl")

include("feature_transform.jl")
include("utils.jl")
using .FeatureTransform

include("deprecations.jl")

export
    # structuring_element.jl
    centered, # from OffsetArrays
    strel,
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
    morpholaplace,

    # connected.jl
    label_components,
    label_components!,
    component_boxes,
    component_lengths,
    component_indices,
    component_subscripts,
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
    clearborder

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
