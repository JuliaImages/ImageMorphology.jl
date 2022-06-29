"""
    SEBox{N}(axes; [r])

The N-dimensional structuring element with all elements connected. This is a special case of
[`SEMask`](@ref) that ImageMorphology algorithms might provide
optimized implementation.

It is recommended to use [`strel_box`](@ref) and [`strel_type`](@ref):

```jldoctest; setup=:(using ImageMorphology; using ImageMorphology.StructuringElements)
julia> se = strel_box((3, 3)) # C8 connectivity
3×3 SEBoxArray{2, UnitRange{$Int}} with indices -1:1×-1:1:
 1  1  1
 1  1  1
 1  1  1

julia> strel_type(se)
SEBox{2, UnitRange{$Int}}((-1:1, -1:1), (1, 1))

julia> se = centered(collect(se)) # converted to normal centered array
3×3 OffsetArray(::Matrix{Bool}, -1:1, -1:1) with eltype Bool with indices -1:1×-1:1:
 1  1  1
 1  1  1
 1  1  1

julia> strel_type(se)
SEMask{2}()
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

"""
    SEBoxArray(se::SEBox)

The instantiated array object of [`SEBox`](@ref SEBox).
"""
struct SEBoxArray{N,R<:AbstractUnitRange{Int}} <: MorphologySEArray{N}
    axes::NTuple{N,R}
    r::Dims{N}
end
SEBoxArray(se::SEBox{N,R}) where {N,R} = SEBoxArray{N,R}(se.axes, se.r)
strel_type(A::SEBoxArray{N}) where {N} = SEBox{N}(A.axes; r=A.r)
strel_size(se::SEBoxArray) = @. 1 + 2 * se.r

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


"""
    strel_box(A; r=1)
    strel_box(size; r=size .÷ 2)

Construct the N-dimensional structuring element (SE) with all elements in the local window
connected.

If image `A` is provided, then the SE size will be `(2r+1, 2r+1, ...)` with default
half-size `r=1`. If `size` is provided, the default `r` will be `size .÷ 2`. The default
`dims` will be all dimensions, that is, `(1, 2, ..., length(size))`.

```jldoctest; setup=:(using ImageMorphology; using ImageMorphology.StructuringElements)
julia> img = rand(64, 64);

julia> strel_box(img)
3×3 SEBoxArray{2, UnitRange{$Int}} with indices -1:1×-1:1:
 1  1  1
 1  1  1
 1  1  1

julia> strel_box(img; r=2)
5×5 SEBoxArray{2, UnitRange{$Int}} with indices -2:2×-2:2:
 1  1  1  1  1
 1  1  1  1  1
 1  1  1  1  1
 1  1  1  1  1
 1  1  1  1  1

julia> strel_box((5,5); r=(1,2))
5×5 SEBoxArray{2, UnitRange{$Int}} with indices -2:2×-2:2:
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
    dims = _to_dims(Val(N), dims)
    sz, r = if isnothing(r)
        ntuple(i -> !isempty(dims) && in(i, dims) ? 3 : 1, N), 1
    elseif r isa Dims{N}
        ntuple(i -> !isempty(dims) && in(i, dims) ? 2r[i] + 1 : 1, N), r
    elseif r isa Integer
        ntuple(i -> !isempty(dims) && in(i, dims) ? 2r + 1 : 1, N), r
    end
    return strel_box(sz, dims)
end
function strel_box(sz::Dims{N}, dims=ntuple(identity, N); r::Union{Nothing,Dims{N},Int}=nothing) where {N}
    dims = _to_dims(Val(N), dims)
    all(isodd, sz) || throw(ArgumentError("size should be odd integers"))
    radius = if isnothing(r)
        ntuple(i -> !isempty(dims) && in(i, dims) ? sz[i] ÷ 2 : 0, N)
    elseif r isa Dims{N}
        r
    elseif r isa Integer
        ntuple(i -> !isempty(dims) && in(i, dims) ? r : 0, N)
    end
    ax = map(r -> (-r):r, sz .÷ 2)
    return SEBoxArray(SEBox{N}(ax; r=radius))
end
