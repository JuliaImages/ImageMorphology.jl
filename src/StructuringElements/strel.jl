"""
    strel([T], X::AbstractArray)

Convert structuring element (SE) `X` to appropriate presentation format with element type `T`.
This is a useful tool to generate SE that most ImageMorphology functions understand.

ImageMorphology currently supports two commonly used representations:

- `T=CartesianIndex`: offsets to its center point. The output type is
  `Vector{CartesianIndex{N}}`.
- `T=Bool`: connectivity mask where `true` indicates connected to its center point. The
  output type is `BitArray{N}`.

```jldoctest; setup=:(using ImageMorphology)
julia> se_mask = centered(Bool[1 1 0; 1 1 0; 0 0 0]) # connectivity mask
3×3 OffsetArray(::Matrix{Bool}, -1:1, -1:1) with eltype Bool with indices -1:1×-1:1:
 1  1  0
 1  1  0
 0  0  0

julia> se_offsets = strel(CartesianIndex, se_mask) # displacement offsets to its center point
3-element Vector{CartesianIndex{2}}:
 CartesianIndex(-1, -1)
 CartesianIndex(0, -1)
 CartesianIndex(-1, 0)

julia> se = strel(Bool, se_offsets)
3×3 OffsetArray(::BitMatrix, -1:1, -1:1) with eltype Bool with indices -1:1×-1:1:
 1  1  0
 1  1  0
 0  0  0
```

See also [`strel_diamond`](@ref) and [`strel_box`](@ref) for SE constructors for two
special cases.
"""
function strel end

strel(se) = strel(strel_type(se), se)

# convenient user interface without exporting MorphologySE
function strel(::Type{ET}, se::AbstractArray) where {ET<:CartesianIndex}
    return strel(SEOffset{strel_ndims(se)}(), se)
end
strel(::Type{ET}, se::AbstractArray) where {ET<:Bool} = strel(SEMask{strel_ndims(se)}(), se)

# conversion between different SE arrays
function strel(SET::MorphologySE, se::T) where {T}
    strel_type(se) == SET || error("unsupported conversion from type $T to $SET")
    return se
end

function strel(::SEMask{N}, offsets::AbstractArray{CartesianIndex{N}}) where {N}
    isempty(offsets) && return centered(trues(ntuple(_ -> 1, N)))
    mn, mx = extrema(offsets)
    r = ntuple(N) do i
        max(abs(mn.I[i]), abs(mx.I[i]))
    end
    sz = @. 2r + 1
    se = centered(falses(sz))
    se[offsets] .= true
    se[zero(eltype(offsets))] = true # always set center point to true
    return centered(BitArray(OffsetArrays.no_offset_view(se)))
end
strel(::SEMask{N}, mask::AbstractArray{Bool,N}) where {N} = mask

function strel(::SEOffset{N}, connectivity::AbstractArray{Bool,N}) where {N}
    all(isodd, size(connectivity)) || error("`connectivity` must be odd-sized")
    ax = axes(connectivity)
    is_symmetric = all(r -> first(r) == -last(r), ax)
    if !is_symmetric && all(first.(axes(connectivity)) .== 1)
        # To keep consistent with the "kernel" concept in ImageFiltering, we require
        # the connectivity mask to be centered as well.
        # This is a perhaps permanent depwarn to throw friendly message to the user
        # if they're used to use, e.g., `trues(3, 3)` as the input.
        msg = "connectivity mask is expected to be a centered bool array"
        hint = "Do you mean `centered(connectivity)`"
        Base.depwarn("$msg. $hint?", :strel)
    elseif !is_symmetric
        throw(ArgumentError("`connectivity` must be symmetric bool array"))
    end
    connectivity = centered(connectivity)
    # always skip center point
    return [i for i in CartesianIndices(connectivity) if connectivity[i] && !iszero(i)]
end

"""
    strel_type(x)

Infer the structuring element type for `x`.

!!! note "developer note"
    This function is used to dispatch special SE types, e.g., [`SEBoxArray`](@ref), to
    optimized implementation of particular morphology filter. In this sense it is required
    for custom SE array types to define this method.
"""
strel_type(se::MorphologySE) = se
strel_type(::AbstractArray{Bool,N}) where {N} = SEMask{N}()
strel_type(::AbstractVector{CartesianIndex{N}}) where {N} = SEOffset{N}()
strel_type(::CartesianIndices{N}) where {N} = SEOffset{N}()
strel_type(::T) where {T} = error("invalid structuring element data type: $T")

"""
    strel_size(x)

Calculate the minimal block size that contains the structuring element. The result
will be a tuple of odd integers.

```jldoctest; setup=:(using ImageMorphology; using ImageMorphology.StructuringElements)
julia> se = strel_diamond((5, 5); r=1)
5×5 SEDiamondArray{2, 2, UnitRange{$Int}, 0} with indices -2:2×-2:2:
 0  0  0  0  0
 0  0  1  0  0
 0  1  1  1  0
 0  0  1  0  0
 0  0  0  0  0

julia> strel_size(se) # is not (5, 5)
(3, 3)

julia> strel(Bool, strel(CartesianIndex, se)) # because it only checks the minimal enclosing block
3×3 OffsetArray(::BitMatrix, -1:1, -1:1) with eltype Bool with indices -1:1×-1:1:
 0  1  0
 1  1  1
 0  1  0

julia> se = [CartesianIndex(1, 1), CartesianIndex(-2, -2)];

julia> strel_size(se) # is not (4, 4)
(5, 5)

julia> strel(Bool, se) # because the connectivity mask has to be odd size
5×5 OffsetArray(::BitMatrix, -2:2, -2:2) with eltype Bool with indices -2:2×-2:2:
 1  0  0  0  0
 0  0  0  0  0
 0  0  1  0  0
 0  0  0  1  0
 0  0  0  0  0

julia> se = strel_diamond((5, 5), (1, ); r=1)
5×5 SEDiamondArray{2, 1, UnitRange{$Int}, 1} with indices -2:2×-2:2:
 0  0  0  0  0
 0  0  1  0  0
 0  0  1  0  0
 0  0  1  0  0
 0  0  0  0  0

julia> strel_size(se)
(3, 1)
```
"""
strel_size(se) = size(strel(Bool, strel(CartesianIndex, se)))

"""
    strel_ndims(x)::Int

Infer the dimension of the structuring element `x`
"""
strel_ndims(se) = strel_ndims(strel_type(se))
strel_ndims(::MorphologySE{N}) where {N} = N
