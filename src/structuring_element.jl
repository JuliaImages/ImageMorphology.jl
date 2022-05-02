abstract type MorphologySE{N} end

"""
    SEMask{N}()

The N-dimensional structuring element in terms of bool array as mask. Typically, `true`
represents the foreground and `false` represents the background.
"""
struct SEMask{N} <: MorphologySE{N} end

"""
    SEOffset{N}()

The N-dimensional structuring element in terms of displacement offset.
"""
struct SEOffset{N} <: MorphologySE{N} end

"""
    SEDiamond{N}([size], [dims]; r=1)

The N-dimensional diamond shape structuring element.
"""
struct SEDiamond{N,K} <: MorphologySE{N}
    size::Dims{N}
    dims::Dims{K}
    r::Int # radius
    function SEDiamond{N,K}(size::Dims{N}, dims::Dims{K}, r) where {N,K}
        all(isodd, size) || throw(ArgumentError("all size should be odd number"))
        _is_unique_tuple(dims) || throw(ArgumentError("dims should be unique"))
        N >= K || throw(ArgumentError("`size` length should be at least $K"))
        all(i->i<=N, dims) || throw(ArgumentError("all `dims` values should be less than or equal to $N"))
        return new{N,K}(size, dims, r)
    end
end
function SEDiamond{N}(sz=ntuple(_ -> 3, N), dims=ntuple(identity, N); r=1) where {N}
    return SEDiamond{N,length(dims)}(sz, dims, r)
end

"""
    SEWindow{N}([size], [dims]; r=1)

The N-dimensional structuring element with all elements connected.
"""
struct SEWindow{N,K} <: MorphologySE{N}
    size::Dims{N}
    dims::Dims{K}
    r::Int # radius
    function SEWindow{N,K}(size::Dims{N}, dims::Dims{K}, r) where {N,K}
        all(isodd, size) || throw(ArgumentError("all size should be odd number"))
        _is_unique_tuple(dims) || throw(ArgumentError("dims should be unique"))
        N >= K || throw(ArgumentError("`size` length should be at least $K"))
        all(i->i<=N, dims) || throw(ArgumentError("all `dims` values should be less than or equal to $N"))
        return new{N,K}(size, dims, r)
    end
end
function SEWindow{N}(sz=ntuple(_ -> 3, N), dims=ntuple(identity, N); r=1) where {N}
    return SEWindow{N,length(dims)}(sz, dims, r)
end

# Helper array type to build a consistant array interface for SEs
struct SEDiamondArray{N,K,S} <: AbstractArray{Bool,N}
    size::Dims{N}
    dims::Dims{K}
    r::Int # radius
    _center::Dims{N}
    _rdims::Dims{S}
end
function SEDiamondArray(se::SEDiamond{N,K}) where {N,K}
    _center = @. (se.size + 1) ÷ 2
    _rdims = _cal_rdims(Val(N), se.dims)
    return SEDiamondArray{N,K,length(_rdims)}(se.size, se.dims, se.r, _center, _rdims)
end

@inline Base.size(A::SEDiamondArray{N}) where {N} = A.size
@inline Base.IndexStyle(::SEDiamondArray) = IndexCartesian()
Base.@propagate_inbounds function Base.getindex(A::SEDiamondArray{N,K}, inds::Int...) where {N,K}
    # for remaining dimensions, check if it is at the center position
    ri = _tuple_getindex(inds, A._rdims)
    rc = _tuple_getindex(A._center, A._rdims)
    ri == rc || return false
    # for masked dimensions, compare if the city-block distance to center is within radius
    mi = _tuple_getindex(inds, A.dims)
    mc = _tuple_getindex(A._center, A.dims)
    r = A.r
    return ifelse(sum(abs, mi .- mc) > r, false, true)
end

struct SEWindowArray{N,K,S} <: AbstractArray{Bool,N}
    dims::Dims{K}
    size::Dims{N}
    r::Int
    _center::Dims{N}
    _rdims::Dims{S}
end
function SEWindowArray(se::SEWindow{N,K}) where {N,K}
    _center = @. (se.size + 1) ÷ 2
    _rdims = _cal_rdims(Val(N), se.dims)
    return SEWindowArray{N,K,length(_rdims)}(se.dims, se.size, se.r, _center, _rdims)
end

@inline Base.size(A::SEWindowArray) = A.size
@inline function Base.getindex(A::SEWindowArray, inds::Int...)
    # for remaining dimensions, check if it is at the center position
    ri = _tuple_getindex(inds, A._rdims)
    rc = _tuple_getindex(A._center, A._rdims)
    ri == rc || return false

    # for masked dimensions, compare if any of the index is within radius
    mi = _tuple_getindex(inds, A.dims)
    mc = _tuple_getindex(A._center, A.dims)
    return ifelse(any(abs.(mi .- mc) .> A.r), false, true)
end

_tuple_getindex(t::Tuple, inds::Dims) = ntuple(i->t[inds[i]], length(inds))
_cal_rdims(::Val{N}, dims::NTuple{K}) where {N,K} = Dims{N-K}(filter(i->!in(i, dims), 1:N))
_is_unique_tuple(t::Tuple) = any(i->t[i] in t[1:i-1], 2:length(t)) ? false : true


"""
    strel_type(x)

Infer the structuring element type for `x`.
"""
strel_type(se::MorphologySE) = se
strel_type(::AbstractArray{Bool,N}) where {N} = SEMask{N}()
strel_type(::AbstractVector{CartesianIndex{N}}) where {N} = SEOffset{N}()
strel_type(::CartesianIndices{N}) where {N} = SEOffset{N}()
strel_type(A::SEDiamondArray{N}) where {N} = SEDiamond{N}(size(A), A.dims; r=A.r)
strel_type(A::SEWindowArray{N}) where {N} = SEWindow{N}(size(A), A.dims)
strel_type(::T) where T = error("invalid structuring element data type: $T")

"""
    strel_size(x)

Calculate the minimal block size that contains the structuring element. The result
will be a tuple of odd integers.

```jldoctest; setup=:(using ImageMorphology)
julia> se = strel_diamond((5, 5))
5×5 ImageMorphology.SEDiamondArray{2, 2, 0}:
 0  0  0  0  0
 0  0  1  0  0
 0  1  1  1  0
 0  0  1  0  0
 0  0  0  0  0

julia> strel_size(se) # is not (5, 5)
(3, 3)

julia> strel(Bool, strel(CartesianIndex, se)) # because it only checks the minimal enclosing block
3×3 BitMatrix:
 0  1  0
 1  1  1
 0  1  0

julia> se = [CartesianIndex(1, 1), CartesianIndex(-2, -2)];

julia> strel_size(se) # is not (4, 4)
(5, 5)

julia> strel(Bool, se) # because the connectivity mask has to be odd size
5×5 BitMatrix:
 1  0  0  0  0
 0  0  0  0  0
 0  0  1  0  0
 0  0  0  1  0
 0  0  0  0  0

julia> se = strel_diamond((5, 5), (1, ))
5×5 ImageMorphology.SEDiamondArray{2, 1, 1}:
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
strel_size(se::SEDiamondArray) = ntuple(i->in(i, se.dims) ? 1+2*se.r : 1, strel_ndims(se))
strel_size(se::SEWindowArray) = ntuple(i->in(i, se.dims) ? 1+2*se.r : 1, strel_ndims(se))

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
julia> se_mask = Bool[1 1 0; 1 1 0; 0 0 0] # connectivity mask
3×3 Matrix{Bool}:
 1  1  0
 1  1  0
 0  0  0

julia> se_offsets = strel(CartesianIndex, se_mask) # displacement offsets to its center point
3-element Vector{CartesianIndex{2}}:
 CartesianIndex(-1, -1)
 CartesianIndex(0, -1)
 CartesianIndex(-1, 0)

julia> se = strel(Bool, se_offsets)
3×3 BitMatrix:
 1  1  0
 1  1  0
 0  0  0
```

See also [`strel_diamond`](@ref) and [`strel_window`](@ref) for SE constructors for two
special cases.
"""
function strel end

strel(se) = strel(strel_type(se), se)

# convenient user interface without exporting MorphologySE
strel(::Type{ET}, se::AbstractArray) where {ET<:CartesianIndex} = strel(SEOffset{strel_ndims(se)}(), se)
strel(::Type{ET}, se::AbstractArray) where {ET<:Bool} = strel(SEMask{strel_ndims(se)}(), se)

# constructor for special SEs
_strel_array(se::SEDiamond) = SEDiamondArray(se)
_strel_array(se::SEWindow) = SEWindowArray(se)

"""
    strel_diamond(img, [dims]; r=1)
    strel_diamond(size, [dims]; r=1)

Construct the N-dimensional structuring element (SE) for a diamond shape. If image is
provided, then `size=(3, 3, ...)` and `(1, 2, ..., N)` are the default values for `size` and
`dims`.

```jldoctest; setup=:(using ImageMorphology)
julia> img = rand(64, 64);

julia> se = strel_diamond(img)
3×3 ImageMorphology.SEDiamondArray{2, 2, 0}:
 0  1  0
 1  1  1
 0  1  0

julia> se = strel_diamond((3,3), (1,)) # 3×3 mask along dimension 1
3×3 ImageMorphology.SEDiamondArray{2, 1, 1}:
 0  1  0
 0  1  0
 0  1  0

julia> se = strel_diamond((3,5); r=2) # 3×5 mask with radius 2
3×5 ImageMorphology.SEDiamondArray{2, 2, 0}:
 0  1  1  1  0
 1  1  1  1  1
 0  1  1  1  0
```

!!! note "specialization and performance"
    The diamond shape `SEDiamond` is a special type that many morpholoy algorithms may
    provide much more efficient implementations for. For this reason, if one tries to
    collect an `SEDiamondArray` into other array types (e.g., `Array{Bool}` via `collect`),
    then a significant performance drop might be very likely to happen.

See also [`strel`](@ref) and [`strel_window`](@ref).
"""
function strel_diamond(img::AbstractArray{T,N}, dims=coords_spatial(img); kw...) where {T,N}
    return strel_diamond(ntuple(i->3, N), dims; kw...)
end
function strel_diamond(sz::Dims{N}, dims::Dims{K}=ntuple(identity, N); kw...) where {N,K}
    return _strel_array(SEDiamond{N}(sz, dims; kw...))
end

"""
    strel_window(img, [dims]; r=1)
    strel_window(size, [dims]; r=1)

Construct the N-dimensional structuring element (SE) with all elements in the local window
connected. If image is provided, then `size=(3, 3, ...)` and `(1, 2, ..., N)` are the
default values for `size` and `dims`.

```jldoctest; setup=:(using ImageMorphology)
julia> img = rand(64, 64);

julia> se = strel_window(img)
3×3 ImageMorphology.SEWindowArray{2, 2, 0}:
 1  1  1
 1  1  1
 1  1  1

julia> se = strel_window((3,3), (1,)) # 3×3 mask along dimension 1
3×3 ImageMorphology.SEWindowArray{2, 1, 1}:
 0  1  0
 0  1  0
 0  1  0

julia> se = strel_window((3,5); r=2) # 3×5 mask with radius 2
3×5 ImageMorphology.SEWindowArray{2, 2, 0}:
 1  1  1  1  1
 1  1  1  1  1
 1  1  1  1  1
```

!!! note "specialization and performance"
    The window shape `SEWindow` is a special type that many morpholoy algorithms may
    provide efficient implementations for. For this reason, if one tries to collect an
    `SEWindowArray` into other array types (e.g., `Array{Bool}` via `collect`), then a
    significant performance drop might be very likely to happen.

See also [`strel`](@ref) and [`strel_window`](@ref).
"""
function strel_window(img::AbstractArray{T,N}, dims=coords_spatial(img); kw...) where {T,N}
    return strel_window(ntuple(i->3, N), dims; kw...)
end
function strel_window(sz::Dims{N}, dims::Dims{K}=ntuple(identity, N); kw...) where {N,K}
    return _strel_array(SEWindow{N}(sz, dims; kw...))
end

# conversion between different SE arrays
function strel(SET::MorphologySE, se::T) where {T}
    strel_type(se) == SET || error("unsupported conversion from type $T to $SET")
    return se
end

function strel(::SEMask{N}, offsets::AbstractArray{CartesianIndex{N}}) where {N}
    isempty(offsets) && return trues(ntuple(_->1, N))
    mn, mx = extrema(offsets)
    r = ntuple(N) do i
        max(abs(mn.I[i]), abs(mx.I[i]))
    end
    sz = @. 2r + 1
    se = OffsetArrays.centered(falses(sz))
    se[offsets] .= true
    se[zero(eltype(offsets))] = true # always set center point to true
    return BitArray(OffsetArrays.no_offset_view(se))
end
strel(::SEMask{N}, mask::AbstractArray{Bool,N}) where {N} = mask

function strel(::SEOffset{N}, connectivity::AbstractArray{Bool,N}) where {N}
    all(isodd, size(connectivity)) || error("`connectivity` must be odd-sized")
    connectivity = OffsetArrays.centered(connectivity)
    # always skip center point
    return [i for i in CartesianIndices(connectivity) if connectivity[i] && !iszero(i)]
end
