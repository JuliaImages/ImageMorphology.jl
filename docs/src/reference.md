# [Reference](@id reference_index)

## Interface overview

The following tables give an overview of ImageMorphology interfaces

[**Structuring Element (SE)**](@ref reference_se)

| name                          | summary |
| :---------------------------- | :------ |
| [`strel`](@ref)               | convert between different SE representations    |
| [`strel_type`](@ref)          | infer the SE type                               |
| [`strel_size`](@ref)          | get the minimal block size that contains the SE |
| [`strel_box`](@ref)           | construct a box-shaped SE, e.g., C8, C26 connectivity |
| [`strel_diamond`](@ref)       | construct a diamond-shaped SE, e.g., C4, C6 connectivity |
| [`OffsetArrays.centered`](@ref OffsetArrays.centered) | shift the array center to `(0, 0, ..., 0)`    |


[**Morphological Operations**](@ref reference_ops)

| name                                                   | summary | examples |
| :----------------------------------------------------- | :------ | ---- |
| [`extreme_filter`](@ref) and [`extreme_filter!`](@ref) | iteratively apply a select function `f(x, y)` to each neighborhood | [`extreme_filter` operation](@ref op_extreme_filter) |
| [`dilate`](@ref) and [`dilate!`](@ref)                 | morphological max filter  | [`dilate` operation](@ref op_dilate)   |
| [`erode`](@ref) and [`erode!`](@ref)                   | morphological min filter  | [`erode` operation](@ref op_erode)     |
| [`opening`](@ref) and [`opening!`](@ref)               | fills white holes         | [`opening` operation](@ref op_opening) |
| [`closing`](@ref) and [`closing!`](@ref)               | fills black holes         | [`closing` operation](@ref op_closing) |
| [`bothat`](@ref) and [`bothat!`](@ref)                 | extract black details     | [`bothat` operation](@ref op_bothat)   |
| [`tophat`](@ref) and [`tophat!`](@ref)                 | extract white details     | [`tophat` operation](@ref op_tophat)   |
| [`morphogradient`](@ref)                               | morphological gradient    |                                        |
| [`morpholaplace`](@ref)                                | morpholigical laplacian   |                                        |

## [Structuring element](@id reference_se)

Structuring element is the key concept in morphology. If you're not familiar with this, please
read [concept: structuring element](@ref concept_se) first.

```@docs
# conversion
strel

# constructor
strel_box
strel_diamond

# helpers
strel_type
strel_size
ImageMorphology.strel_ndims
OffsetArrays.centered
OffsetArrays.center

ImageMorphology.is_symmetric

## types
ImageMorphology.SEMask
ImageMorphology.SEOffset
ImageMorphology.SEDiamond
ImageMorphology.SEBox
ImageMorphology.SEDiamondArray
ImageMorphology.SEBoxArray
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
morphogradient
morpholaplace
```

## Components and segmentations

```@docs
label_components
component_boxes
component_lengths
component_indices
component_subscripts
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
