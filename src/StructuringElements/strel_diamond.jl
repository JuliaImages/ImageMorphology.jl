"""
    SEDiamond{N}(axes, [dims]; [r])

A (holy) trait type for the N-dimensional diamond shape structuring element. This is a
special case of [`SEMask`](@ref ImageMorphology.SEMask) that ImageMorphology algorithms
might provide optimized implementation.

It is recommended to use [`strel_diamond`](@ref) and [`strel_type`](@ref):

```jldoctest; setup=:(using ImageMorphology; using ImageMorphology.StructuringElements)
julia> se = strel_diamond((3, 3)) # C4 connectivity
3×3 SEDiamondArray{2, 2, UnitRange{$Int}, 0} with indices -1:1×-1:1:
 0  1  0
 1  1  1
 0  1  0

julia> strel_type(se)
SEDiamond{2, 2, UnitRange{$Int}}((-1:1, -1:1), (1, 2), 1)

julia> se = centered(collect(se)) # converted to normal centered array
3×3 OffsetArray(::Matrix{Bool}, -1:1, -1:1) with eltype Bool with indices -1:1×-1:1:
 0  1  0
 1  1  1
 0  1  0

julia> strel_type(se)
SEMask{2}()
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
strel_type(A::SEDiamondArray{N}) where {N} = SEDiamond{N}(A.axes, A.dims; r=A.r)
function strel_size(se::SEDiamondArray)
    return ntuple(i -> in(i, se.dims) ? 1 + 2 * se.r : 1, strel_ndims(se))
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

_tuple_getindex(t::Tuple, inds::Dims) = ntuple(i -> t[inds[i]], length(inds))
function _cal_rdims(::Val{N}, dims::NTuple{K}) where {N,K}
    return Dims{N - K}(filter(i -> !in(i, dims), 1:N))
end
_is_unique_tuple(t::Tuple) = any(i -> t[i] in t[1:(i - 1)], 2:length(t)) ? false : true

"""
    strel_diamond(A::AbstractArray, [dims]; r=1)
    strel_diamond(size, [dims]; [r])

Construct the N-dimensional structuring element (SE) for a diamond shape.

If image `A` is provided, then the SE size will be `(2r+1, 2r+1, ...)` with default
half-size `r=1`. If `size` is provided, the default `r` will be `maximum(size)÷2`. The
default `dims` will be all dimensions, that is, `(1, 2, ..., length(size))`.

```jldoctest; setup=:(using ImageMorphology; using ImageMorphology.StructuringElements)
julia> img = rand(64, 64);

julia> strel_diamond(img) # default size for image input is (3, 3)
3×3 SEDiamondArray{2, 2, UnitRange{$Int}, 0} with indices -1:1×-1:1:
 0  1  0
 1  1  1
 0  1  0

julia> strel_diamond(img; r=2) # equivalent to `strel_diamond((5,5))`
5×5 SEDiamondArray{2, 2, UnitRange{$Int}, 0} with indices -2:2×-2:2:
 0  0  1  0  0
 0  1  1  1  0
 1  1  1  1  1
 0  1  1  1  0
 0  0  1  0  0

julia> strel_diamond(img, (1,)) # mask along dimension 1
3×1 SEDiamondArray{2, 1, UnitRange{$Int}, 1} with indices -1:1×0:0:
 1
 1
 1

julia> strel_diamond((3,3), (1,)) # 3×3 mask along dimension 1
3×3 SEDiamondArray{2, 1, UnitRange{$Int}, 1} with indices -1:1×-1:1:
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
    return SEDiamondArray(SEDiamond{N}(ax, dims; kw...))
end

# Tuple(1) is not inferable
@inline _to_dims(i::Int) = (i,)
@inline _to_dims(dims::Dims) = dims
@inline _to_dims(v) = Tuple(v) # fallback
