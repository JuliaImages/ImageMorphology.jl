const _docstring_se = """
`se` is the structuring element that defines the neighborhood of the image. See
[`strel`](@ref) for more details. If `se` is not specified, then it will use the
[`strel_box`](@ref) with an extra keyword `dims` to control the dimensions to filter,
and half-size `r` to control the diamond size.
"""

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
ImageMorphology.SEMask{2}()
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
ImageMorphology.SEOffset{2}()
```

See also [`SEMask`](@ref ImageMorphology.SEMask) for the connectivity mask representation.
More details can be found on he documentation page [Structuring Element](@ref concept_se).
"""
struct SEOffset{N} <: MorphologySE{N} end

"""
    SEDiamond{N}(axes, [dims]; [r])

A (holy) trait type for the N-dimensional diamond shape structuring element. This is a
special case of [`SEMask`](@ref ImageMorphology.SEMask) that ImageMorphology algorithms
might provide optimized implementation.

It is recommended to use [`strel_diamond`](@ref) and [`strel_type`](@ref):

```jldoctest; setup=:(using ImageMorphology)
julia> using OffsetArrays: centered

julia> se = strel_diamond((3, 3)) # C4 connectivity
3×3 ImageMorphology.SEDiamondArray{2, 2, UnitRange{$Int}, 0} with indices -1:1×-1:1:
 0  1  0
 1  1  1
 0  1  0

julia> strel_type(se)
ImageMorphology.SEDiamond{2, 2, UnitRange{$Int}}((-1:1, -1:1), (1, 2), 1)

julia> se = centered(collect(se)) # converted to normal centered array
3×3 OffsetArray(::Matrix{Bool}, -1:1, -1:1) with eltype Bool with indices -1:1×-1:1:
 0  1  0
 1  1  1
 0  1  0

julia> strel_type(se)
ImageMorphology.SEMask{2}()
```
"""
struct SEDiamond{N,K,R<:AbstractUnitRange{Int}} <: MorphologySE{N}
    axes::NTuple{N,R}
    dims::Dims{K}
    r::Int # radius
    function SEDiamond{N,K,R}(axes::NTuple{N,R}, dims::Dims{K}, r) where {N,K,R<:AbstractUnitRange{Int}}
        if !all(r -> first(r) == -last(r), axes)
            throw(ArgumentError("axes must be symmetric along each dimension"))
        end
        _is_unique_tuple(dims) || throw(ArgumentError("dims should be unique"))
        N >= K || throw(ArgumentError("`axes` length should be at least $K"))
        if !all(i -> i <= N, dims)
            throw(ArgumentError("all `dims` values should be less than or equal to $N"))
        end
        return new{N,K,R}(axes, dims, r)
    end
end
function SEDiamond{N}(ax::NTuple{N,R}, dims=ntuple(identity, N); r=maximum(length.(ax)) ÷ 2) where {N,R}
    return SEDiamond{N,length(dims),R}(ax, dims, r)
end

"""
    SEBox{N}(axes; [r])

The N-dimensional structuring element with all elements connected. This is a special case of
[`SEMask`](@ref ImageMorphology.SEMask) that ImageMorphology algorithms might provide
optimized implementation.

It is recommended to use [`strel_box`](@ref) and [`strel_type`](@ref):

```jldoctest; setup=:(using ImageMorphology)
julia> using OffsetArrays: centered

julia> se = strel_box((3, 3)) # C8 connectivity
3×3 ImageMorphology.SEBoxArray{2, UnitRange{$Int}} with indices -1:1×-1:1:
 1  1  1
 1  1  1
 1  1  1

julia> strel_type(se)
ImageMorphology.SEBox{2, UnitRange{$Int}}((-1:1, -1:1), (1, 1))

julia> se = centered(collect(se)) # converted to normal centered array
3×3 OffsetArray(::Matrix{Bool}, -1:1, -1:1) with eltype Bool with indices -1:1×-1:1:
 1  1  1
 1  1  1
 1  1  1

julia> strel_type(se)
ImageMorphology.SEMask{2}()
```
"""
struct SEBox{N,R} <: MorphologySE{N}
    axes::NTuple{N,R}
    r::Dims{N}
    function SEBox{N,R}(axes::NTuple{N,R}, r::Dims{N}) where {N,R<:AbstractUnitRange{Int}}
        if !all(r -> first(r) == -last(r), axes)
            throw(ArgumentError("axes must be symmetric along each dimension"))
        end
        return new{N,R}(axes, r)
    end
end
function SEBox{N}(ax::NTuple{N,R}; r=map(R -> length(R) ÷ 2, ax)) where {N,R}
    r = r isa Integer ? ntuple(_ -> r, N) : r
    return SEBox{N,R}(ax, r)
end

# Helper array type to build a consistant array interface for SEs
"""
    SEDiamondArray(se::SEDiamond)

The instantiated array object of [`SEDiamond`](@ref ImageMorphology.SEDiamond).
"""
struct SEDiamondArray{N,K,R<:AbstractUnitRange{Int},S} <: MorphologySEArray{N}
    axes::NTuple{N,R}
    dims::Dims{K}
    r::Int # radius
    _rdims::Dims{S}
end
function SEDiamondArray(se::SEDiamond{N,K,R}) where {N,K,R}
    _rdims = _cal_rdims(Val(N), se.dims)
    return SEDiamondArray{N,K,R,length(_rdims)}(se.axes, se.dims, se.r, _rdims)
end

@inline Base.axes(A::SEDiamondArray) = A.axes
@inline Base.size(A::SEDiamondArray) = map(length, axes(A))
@inline Base.IndexStyle(::SEDiamondArray) = IndexCartesian()
Base.@propagate_inbounds function Base.getindex(A::SEDiamondArray{N,K}, inds::Int...) where {N,K}
    # for remaining dimensions, check if it is at the center position
    ri = _tuple_getindex(inds, A._rdims)
    all(iszero, ri) || return false
    # for masked dimensions, compare if the city-block distance to center is within radius
    mi = _tuple_getindex(inds, A.dims)
    return ifelse(sum(abs, mi) > A.r, false, true)
end

"""
    SEBoxArray(se::SEBox)

The instantiated array object of [`SEBox`](@ref ImageMorphology.SEBox).
"""
struct SEBoxArray{N,R<:AbstractUnitRange{Int}} <: MorphologySEArray{N}
    axes::NTuple{N,R}
    r::Dims{N}
end
SEBoxArray(se::SEBox{N,R}) where {N,R} = SEBoxArray{N,R}(se.axes, se.r)

@inline Base.axes(A::SEBoxArray) = A.axes
@inline Base.size(A::SEBoxArray) = map(length, axes(A))
@inline function Base.getindex(A::SEBoxArray, inds::Int...)
    @inbounds for i in 1:length(inds)
        if abs(inds[i]) > A.r[i]
            return false
        end
    end
    return true
end

_tuple_getindex(t::Tuple, inds::Dims) = ntuple(i -> t[inds[i]], length(inds))
function _cal_rdims(::Val{N}, dims::NTuple{K}) where {N,K}
    return Dims{N - K}(filter(i -> !in(i, dims), 1:N))
end
_is_unique_tuple(t::Tuple) = any(i -> t[i] in t[1:(i - 1)], 2:length(t)) ? false : true


"""
    strel_type(x)

Infer the structuring element type for `x`.
"""
strel_type(se::MorphologySE) = se
strel_type(::AbstractArray{Bool,N}) where {N} = SEMask{N}()
strel_type(::AbstractVector{CartesianIndex{N}}) where {N} = SEOffset{N}()
strel_type(::CartesianIndices{N}) where {N} = SEOffset{N}()
strel_type(A::SEDiamondArray{N}) where {N} = SEDiamond{N}(A.axes, A.dims; r=A.r)
strel_type(A::SEBoxArray{N}) where {N} = SEBox{N}(A.axes; r=A.r)
strel_type(::T) where {T} = error("invalid structuring element data type: $T")

"""
    strel_size(x)

Calculate the minimal block size that contains the structuring element. The result
will be a tuple of odd integers.

```jldoctest; setup=:(using ImageMorphology)
julia> se = strel_diamond((5, 5); r=1)
5×5 ImageMorphology.SEDiamondArray{2, 2, UnitRange{$Int}, 0} with indices -2:2×-2:2:
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
5×5 ImageMorphology.SEDiamondArray{2, 1, UnitRange{$Int}, 1} with indices -2:2×-2:2:
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
function strel_size(se::SEDiamondArray)
    return ntuple(i -> in(i, se.dims) ? 1 + 2 * se.r : 1, strel_ndims(se))
end
strel_size(se::SEBoxArray) = @. 1 + 2 * se.r

"""
    strel_ndims(x)::Int

Infer the dimension of the structuring element `x`
"""
strel_ndims(se) = strel_ndims(strel_type(se))
strel_ndims(::MorphologySE{N}) where {N} = N

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

# constructor for special SEs
_strel_array(se::SEDiamond) = SEDiamondArray(se)
_strel_array(se::SEBox) = SEBoxArray(se)

"""
    strel_diamond(A::AbstractArray, [dims]; r=1)
    strel_diamond(size, [dims]; [r])

Construct the N-dimensional structuring element (SE) for a diamond shape.

If image `A` is provided, then the SE size will be `(2r+1, 2r+1, ...)` with default
half-size `r=1`. If `size` is provided, the default `r` will be `maximum(size)÷2`. The
default `dims` will be all dimensions, that is, `(1, 2, ..., length(size))`.

```jldoctest; setup=:(using ImageMorphology)
julia> img = rand(64, 64);

julia> strel_diamond(img) # default size for image input is (3, 3)
3×3 ImageMorphology.SEDiamondArray{2, 2, UnitRange{$Int}, 0} with indices -1:1×-1:1:
 0  1  0
 1  1  1
 0  1  0

julia> strel_diamond(img; r=2) # equivalent to `strel_diamond((5,5))`
5×5 ImageMorphology.SEDiamondArray{2, 2, UnitRange{$Int}, 0} with indices -2:2×-2:2:
 0  0  1  0  0
 0  1  1  1  0
 1  1  1  1  1
 0  1  1  1  0
 0  0  1  0  0

julia> strel_diamond(img, (1,)) # mask along dimension 1
3×1 ImageMorphology.SEDiamondArray{2, 1, UnitRange{$Int}, 1} with indices -1:1×0:0:
 1
 1
 1

julia> strel_diamond((3,3), (1,)) # 3×3 mask along dimension 1
3×3 ImageMorphology.SEDiamondArray{2, 1, UnitRange{$Int}, 1} with indices -1:1×-1:1:
 0  1  0
 0  1  0
 0  1  0
```

!!! note "specialization and performance"
    The diamond shape `SEDiamond` is a special type for which many morphology algorithms may
    provide much more efficient implementations. For this reason, if one tries to collect an
    `SEDiamondArray` into other array types (e.g. `Array{Bool}` via `collect`), then a
    significant performance drop is very likely to occur.

See also [`strel`](@ref) and [`strel_box`](@ref).
"""
function strel_diamond(img::AbstractArray{T,N}, dims=coords_spatial(img); r::Union{Nothing,Int}=nothing) where {T,N}
    dims = _to_dims(dims)
    sz, r = if isnothing(r)
        ntuple(i -> in(i, dims) ? 3 : 1, N), 1
    else
        ntuple(i -> in(i, dims) ? 2r + 1 : 1, N), r
    end
    return strel_diamond(sz, dims; r)
end
function strel_diamond(sz::Dims{N}, dims=ntuple(identity, N); kw...) where {N}
    dims = _to_dims(dims)
    all(isodd, sz) || throw(ArgumentError("size should be odd integers"))
    ax = map(r -> (-r):r, sz .÷ 2)
    return _strel_array(SEDiamond{N}(ax, dims; kw...))
end

"""
    strel_box(A; r=1)
    strel_box(size; r=size .÷ 2)

Construct the N-dimensional structuring element (SE) with all elements in the local window
connected.

If image `A` is provided, then the SE size will be `(2r+1, 2r+1, ...)` with default
half-size `r=1`. If `size` is provided, the default `r` will be `size .÷ 2`. The default
`dims` will be all dimensions, that is, `(1, 2, ..., length(size))`.

```jldoctest; setup=:(using ImageMorphology)
julia> img = rand(64, 64);

julia> strel_box(img)
3×3 ImageMorphology.SEBoxArray{2, UnitRange{$Int}} with indices -1:1×-1:1:
 1  1  1
 1  1  1
 1  1  1

julia> strel_box(img; r=2)
5×5 ImageMorphology.SEBoxArray{2, UnitRange{$Int}} with indices -2:2×-2:2:
 1  1  1  1  1
 1  1  1  1  1
 1  1  1  1  1
 1  1  1  1  1
 1  1  1  1  1

julia> strel_box((5,5); r=(1,2))
5×5 ImageMorphology.SEBoxArray{2, UnitRange{$Int}} with indices -2:2×-2:2:
 0  0  0  0  0
 1  1  1  1  1
 1  1  1  1  1
 1  1  1  1  1
 0  0  0  0  0
```

!!! note "specialization and performance"
    The box shape `SEBox` is a special type for which many morphology algorithms may provide
    efficient implementations. For this reason, if one tries to collect an `SEBoxArray` into
    other array types (e.g. `Array{Bool}` via `collect`), then a significant performance
    drop is very likely to occur.

See also [`strel`](@ref) and [`strel_box`](@ref).
"""
function strel_box(A::AbstractArray{T,N}, dims=coords_spatial(A); r::Union{Nothing,Dims{N},Int}=nothing) where {T,N}
    sz, r = if isnothing(r)
        ntuple(i -> in(i, dims) ? 3 : 1, N), 1
    elseif r isa Dims{N}
        ntuple(i -> in(i, dims) ? 2r[i] + 1 : 1, N), r
    elseif r isa Integer
        ntuple(i -> in(i, dims) ? 2r + 1 : 1, N), r
    end
    return strel_box(sz, dims)
end
function strel_box(sz::Dims{N}, dims=ntuple(identity, N); r::Union{Nothing,Dims{N},Int}=nothing) where {N}
    dims = _to_dims(dims)
    all(isodd, sz) || throw(ArgumentError("size should be odd integers"))
    radius = if isnothing(r)
        ntuple(i -> in(i, dims) ? sz[i] ÷ 2 : 0, N)
    elseif r isa Dims{N}
        r
    elseif r isa Integer
        ntuple(i -> in(i, dims) ? r : 0, N)
    end
    ax = map(r -> (-r):r, sz .÷ 2)
    return _strel_array(SEBox{N}(ax; r=radius))
end

# Tuple(1) is not inferable
@inline _to_dims(i::Int) = (i,)
@inline _to_dims(dims::Dims) = dims
@inline _to_dims(v) = Tuple(v) # fallback

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
