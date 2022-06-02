module StructuringElements

export strel, strel_type, strel_size, strel_ndims
export SEMask, MorphologySE, SEOffset
export strel_box, SEBox, SEBoxArray
export strel_diamond, SEDiamond, SEDiamondArray

# TODO(johnnychen94): remove ImageCore dependency
using ImageCore: coords_spatial

using OffsetArrays
using OffsetArrays: centered

abstract type MorphologySE{N} end
abstract type MorphologySEArray{N} <: AbstractArray{Bool,N} end

OffsetArrays.centered(A::MorphologySEArray) = A

"""
    SEMask{N}()

A (holy) trait type for representing structuring element as connectivity mask. This
connectivity mask SE is a bool array where `true` indicates that pixel position is connected
to the center point.

```jldoctest; setup=:(using ImageMorphology)
julia> se = centered(Bool[0 1 0; 1 1 1; 0 1 0]) # commonly known as C4 connectivity
3×3 OffsetArray(::Matrix{Bool}, -1:1, -1:1) with eltype Bool with indices -1:1×-1:1:
 0  1  0
 1  1  1
 0  1  0

julia> strel_type(se)
ImageMorphology.StructuringElements.SEMask{2}()
```

See also [`SEOffset`](@ref ImageMorphology.SEOffset) for the displacement offset
representation. More details can be found on he documentation page [Structuring
Element](@ref concept_se).
"""
struct SEMask{N} <: MorphologySE{N} end

"""
    SEOffset{N}()

A (holy) trait type for representing structuring element as displacement offsets. This
displacement offsets SE is an array of `CartesianIndex` where each element stores the
displacement offset from the center point.

```jldoctest; setup=:(using ImageMorphology)
julia> se = [CartesianIndex(-1, 0), CartesianIndex(0, -1), CartesianIndex(1, 0), CartesianIndex(0, 1)]
4-element Vector{CartesianIndex{2}}:
 CartesianIndex(-1, 0)
 CartesianIndex(0, -1)
 CartesianIndex(1, 0)
 CartesianIndex(0, 1)

julia> strel_type(se)
ImageMorphology.StructuringElements.SEOffset{2}()
```

See also [`SEMask`](@ref ImageMorphology.SEMask) for the connectivity mask representation.
More details can be found on he documentation page [Structuring Element](@ref concept_se).
"""
struct SEOffset{N} <: MorphologySE{N} end

include("strel.jl")
# SE constructors
include("strel_box.jl")
include("strel_diamond.jl")

end
