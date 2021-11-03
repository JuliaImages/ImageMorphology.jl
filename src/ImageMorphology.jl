module ImageMorphology

using ImageCore
using ImageCore: GenericGrayImage
using LinearAlgebra
using TiledIteration: EdgeIterator, SplitAxis, SplitAxes
using Requires

include("convexhull.jl")
include("connected.jl")
include("clearborder.jl")
include("dilation_and_erosion.jl")
include("find_boundaries.jl")
include("thinning.jl")
include("imfill.jl")
include("maxtree.jl")

include("feature_transform.jl")
using .FeatureTransform

include("deprecations.jl")

export
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

    # find_boundaries.jl
    find_boundaries,
    find_boundaries!,

    # thinning.jl
    thinning,
    GuoAlgo,
    imfill,

    # maxtree.jl
    MaxTree,
    areas, boundingboxes, diameters,
    area_opening, area_opening!, area_closing, area_closing!,
    diameter_opening, diameter_opening!, diameter_closing, diameter_closing!,
    local_maxima!, local_maxima, local_minima!, local_minima,

    #feature_transform.jl
    feature_transform,
    distance_transform,

    clearborder

function __init__()
    @require ImageMetadata = "bc367c6b-8a6b-528e-b4bd-a4b897500b49" begin
        # morphological operations for ImageMeta
        dilate(img::ImageMetadata.ImageMeta, region=coords_spatial(img)) = ImageMetadata.shareproperties(img, dilate!(copy(ImageMetadata.arraydata(img)), region))
        erode(img::ImageMetadata.ImageMeta, region=coords_spatial(img)) = ImageMetadata.shareproperties(img, erode!(copy(ImageMetadata.arraydata(img)), region))
    end
end

end # module
