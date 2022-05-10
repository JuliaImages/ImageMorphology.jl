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
include("dilation_and_erosion.jl")
include("isboundary.jl")
include("thinning.jl")
include("imfill.jl")
include("maxtree.jl")

include("feature_transform.jl")
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
    opening,
    closing,
    tophat,
    bothat,
    morphogradient,
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
            return ImageMetadata.shareproperties(img, dilate!(copy(ImageMetadata.arraydata(img)); kwargs...))
        end
        function erode(img::ImageMetadata.ImageMeta; kwargs...)
            return ImageMetadata.shareproperties(img, erode!(copy(ImageMetadata.arraydata(img)); kwargs...))
        end
    end
end

end # module
