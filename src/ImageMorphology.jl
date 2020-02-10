module ImageMorphology

using ImageCore
using ImageCore: GenericGrayImage
using LinearAlgebra
using ColorVectorSpace
using Compat # for CartesianIndices ranges and oneunit()
using Base.Cartesian # TODO: delete this

include("convexhull.jl")
include("connected.jl")

include("dilation_and_erosion.jl")
include("thinning.jl")
include("imfill.jl")
include("maxtree.jl")

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

    # thinning.jl
    thinning,
    GuoAlgo,
    imfill,

    # maxtree.jl
    MaxTree,
    areas, boundingboxes, diameters,
    area_opening, area_opening!, area_closing, area_closing!,
    diameter_opening, diameter_opening!, diameter_closing, diameter_closing!,
    local_maxima!, local_maxima, local_minima!, local_minima

end # module
