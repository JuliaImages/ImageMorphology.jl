# Reference

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

## types
ImageMorphology.SEMask
ImageMorphology.SEOffset
ImageMorphology.SEDiamond
ImageMorphology.SEBox
ImageMorphology.SEDiamondArray
ImageMorphology.SEBoxArray
```
