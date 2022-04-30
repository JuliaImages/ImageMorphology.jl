"""
    extreme_filter(f, A, [dims]) -> out
    extreme_filter(f, A, Ω) -> out

Filter the array `A` using select function `f(x, y)` for each Ω-neighborhood. The name
"extreme" comes from the fact that typical select function `f` choice is `min` and `max`.

For each pixel `p` in `A`, the select function `f` is applied to its Ω-neighborhood
iteratively in a `f(...(f(f(A[p], A[p+Ω[1]]), A[p+Ω[2]]), ...)` manner. For instance, in the
1-dimensional case, `out[p] = f(f(A[p], A[p-1]), A[p+1])` for each `p` is the default
behavior.

The Ω-neighborhood is defined by the `dims` or `Ω` argument. The `dims` argument specifies
the dimensions of the neighborhood that can be used to construct a diamond shape `Ω` using
[`strel_diamond`](@ref). The `Ω` is also known as structuring element (SE), it can be either
displacement offsets or bool array mask, please refer to [`strel`](@ref) for more details.

```jldoctest; setup=:(using ImageMorphology)
julia> M = [4 6 5 3 4; 8 6 9 4 8; 7 8 4 9 6; 6 2 2 1 7; 1 6 5 2 6]
5×5 Matrix{$Int}:
 4  6  5  3  4
 8  6  9  4  8
 7  8  4  9  6
 6  2  2  1  7
 1  6  5  2  6

julia> extreme_filter(max, M) # max-filter using 4 direct neighbors along both dimensions
5×5 Matrix{$Int}:
 8  6  9  5  8
 8  9  9  9  8
 8  8  9  9  9
 7  8  5  9  7
 6  6  6  6  7

julia> extreme_filter(max, M, (1, )) # max-filter along the first dimension (column)
5×5 Matrix{$Int}:
 8  6  9  4  8
 8  8  9  9  8
 8  8  9  9  8
 7  8  5  9  7
 6  6  5  2  7
```

`Ω` can be either an `AbstractArray{Bool}` mask array with `true` element indicating
connectivity, or a `AbstractArray{<:CartesianIndex}` array with each element indicating the
displacement offset to its center element.

```
julia> Ω_mask = Bool[1 1 0; 1 1 0; 1 0 0] # custom neighborhood in mask format
3×3 Matrix{Bool}:
 1  1  0
 1  1  0
 1  0  0

julia> out = extreme_filter(max, M, Ω_mask)
5×5 Matrix{$Int}:
 4  8  6  9  4
 8  8  9  9  9
 8  8  9  9  9
 7  8  8  9  9
 6  6  6  5  7

julia> Ω_offsets = strel(CartesianIndex, Ω_mask) # custom neighborhood as displacement offset
4-element Vector{CartesianIndex{2}}:
 CartesianIndex(-1, -1)
 CartesianIndex(0, -1)
 CartesianIndex(1, -1)
 CartesianIndex(-1, 0)

julia> out == extreme_filter(max, M, Ω_offsets) # both versions work equivalently
true
```

See also the in-place version [`extreme_filter!`](@ref). Another function in ImageFiltering
package `ImageFiltering.mapwindow` provides similar functionality.
"""
extreme_filter(f, A, dims::Dims) = extreme_filter(f, A, strel_diamond(A, dims))
extreme_filter(f, A, Ω::AbstractArray=strel_diamond(A)) = extreme_filter!(f, similar(A), A, Ω)

"""
    extreme_filter!(f, out, A, [dims])
    extreme_filter!(f, out, A, Ω)

The in-place version of [`extreme_filter`](@ref) where `out` is the output array that gets
modified.
"""
extreme_filter!(f, out, A, dims::Dims) = extreme_filter!(f, out, A, strel_diamond(A, dims))
function extreme_filter!(f, out, A, Ω::AbstractArray=strel_diamond(A))
    _extreme_filter!(strel_type(Ω), f, out, A, Ω)
end

function _extreme_filter!(::MorphologySE, f, out, A, Ω)
    axes(out) == axes(A) || throw(DimensionMismatch("axes(out) must match axes(A)"))

    Ω = strel(CartesianIndex, Ω)
    δ = CartesianIndex(strel_size(Ω) .÷ 2)

    R = CartesianIndices(A)
    R_inner = (first(R)+δ):(last(R)-δ)

    # for interior points, boundary check is unnecessary
    @inbounds for p in R_inner
        out[p] = _select_region(f, A, p, Ω)
    end
    # for edge points, skip if the offset exceeds the boundary
    @inbounds for p in EdgeIterator(R, R_inner)
        Ωp = filter(o->in(p+o, R), Ω)
        out[p] = _select_region(f, A, p, Ωp)
    end
    return out
end

# TODO(johnnychen94): add specialization on max/min for Bool array
function _select_region(f, A, p, offsets)
    # Carefully building `offsets` in the caller `_extreme_filter!` so that boundscheck can
    # be skipped here
    s = @inbounds A[p]
    @inbounds for o in offsets
        v = A[p+o]
        s = f(s, v)
    end
    return s
end
