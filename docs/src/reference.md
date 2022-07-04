# [Reference](@id reference_index)

## [Structuring element](@id reference_se)

```@docs
# conversion
strel

# constructor
strel_chain
strel_product
strel_box
strel_diamond

# helpers
strel_type
strel_size
ImageMorphology.StructuringElements.strel_ndims
ImageMorphology.StructuringElements.strel_split
OffsetArrays.centered
OffsetArrays.center

ImageMorphology.StructuringElements.is_symmetric

## types
ImageMorphology.StructuringElements.SEMask
ImageMorphology.StructuringElements.SEOffset
ImageMorphology.StructuringElements.SEDiamond
ImageMorphology.StructuringElements.SEBox
ImageMorphology.StructuringElements.SEDiamondArray
ImageMorphology.StructuringElements.SEBoxArray
```

## [Morphological operations](@id reference_ops)

```@docs
extreme_filter
extreme_filter!
dilate
dilate!
erode
erode!
opening
opening!
closing
closing!
tophat
tophat!
bothat
bothat!
mgradient
mgradient!
mlaplacian
mlaplacian!
```

## Geodesic operations

```@docs
mreconstruct
mreconstruct!
underbuild
underbuild!
overbuild
overbuild!
```

## Components and segmentation

```@docs
label_components
label_components!
component_boxes
component_lengths
component_indices
component_centroids
```

## Max tree

```@docs
MaxTree
areas
boundingboxes
diameters
area_opening
area_opening!
area_closing
area_closing!
diameter_opening
diameter_opening!
diameter_closing
diameter_closing!
local_maxima!
local_maxima
local_minima!
local_minima
ImageMorphology.rebuild!
ImageMorphology.filter_components!
```

## Feature transform

```@docs
feature_transform
distance_transform
clearborder
```

## Misc

```@docs
convexhull
isboundary
isboundary!
ImageMorphology.isboundary_thick
```
