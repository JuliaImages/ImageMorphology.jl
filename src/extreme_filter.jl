"""
    extreme_filter(f, A; [dims]) -> out
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

# Examples

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

julia> extreme_filter(max, M; dims=(1, )) # max-filter along the first dimension (column)
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
extreme_filter(f, A; dims::Dims=coords_spatial(A)) = extreme_filter(f, A, strel_diamond(A, dims))
extreme_filter(f, A, Ω::AbstractArray) = extreme_filter!(f, similar(A), A, Ω)

"""
    extreme_filter!(f, out, A, [dims])
    extreme_filter!(f, out, A, Ω)

The in-place version of [`extreme_filter`](@ref) where `out` is the output array that gets
modified.
"""
extreme_filter!(f, out, A, dims::Dims) = extreme_filter!(f, out, A, strel_diamond(A, dims))
function extreme_filter!(f, out, A, Ω::AbstractArray=strel_diamond(A))
    axes(out) == axes(A) || throw(DimensionMismatch("axes(out) must match axes(A)"))
    _extreme_filter!(strel_type(Ω), f, out, A, Ω)
end

_extreme_filter!(::MorphologySE, f, out, A, Ω) = _extreme_filter_generic!(f, out, A, Ω)
_extreme_filter!(::SEDiamond, f, out, A, Ω) = _extreme_filter_diamond!(f, out, A, Ω)



###
# Implementation details
###

function _extreme_filter_generic!(f, out, A, Ω)
    @debug "call the generic `extreme_filter` implementation" fname=_extreme_filter_generic!
    Ω = strel(CartesianIndex, Ω)
    δ = CartesianIndex(strel_size(Ω) .÷ 2)

    R = CartesianIndices(A)
    R_inner = (first(R)+δ):(last(R)-δ)

    # NOTE(johnnychen94): delibrately duplicate the kernel codes here for inner loop and
    # boundary loop, because keeping the big function body allows Julia to do further
    # optimization and thus significantly improves the overall performance.
    @inbounds for p in R_inner
        # for interior points, boundary check is unnecessary
        s = A[p]
        for o in Ω
            v = A[p+o]
            s = f(s, v)
        end
        out[p] = s
    end
    @inbounds for p in EdgeIterator(R, R_inner)
        # for edge points, skip if the offset exceeds the boundary
        s = A[p]
        for o in Ω
            q = p+o
            checkbounds(Bool, A, q) || continue
            s = f(s, A[q])
        end
        out[p] = s
    end
    return out
end

# optimized implementation for SEDiamond -- a typical case of separable filter
function _extreme_filter_diamond!(f, out, A, Ω::SEDiamondArray)
    @debug "call the optimized `extreme_filter` implementation for SEDiamond SE" fname=_extreme_filter_diamond!
    rΩ = strel_size(Ω) .÷ 2
    r = maximum(rΩ)

    # To avoid the result affected by loop order, we need two arrays
    src = (out === A) || (r>1) ? copy(A) : A
    out .= src

    # applying radius=r filter is equivalent to applying radius=1 filter r times
    for i in 1:r
        Ω = strel_diamond(A, Ω.dims)
        rΩ = strel_size(Ω) .÷ 2
        inds = axes(A)
        for d in 1:ndims(A)
            # separately apply to each dimension
            if size(out, d) == 1 || rΩ[d] == 0
                continue
            end
            Rpre = CartesianIndices(inds[1:d-1])
            Rpost = CartesianIndices(inds[d+1:end])
            _extreme_filter_C2!(f, out, src, Rpre, inds[d], Rpost)
        end

        if r > 1 && i < r
            src .= out
        end
    end
    return out
end

# fast version for C2 connectivity along each dimension
@noinline function _extreme_filter_C2!(f, dst, src, Rpre, inds, Rpost)
    # manually unroll the loop for r=1 case
    # for r>1 case, it is equivalently to apply r times
    @inbounds for Ipost in Rpost, Ipre in Rpre
        # loop inner region
        for i in first(inds)+1:last(inds)-1
            a1 = src[Ipre, i-1, Ipost]
            a2 = dst[Ipre, i, Ipost]
            a3 = src[Ipre, i+1, Ipost]
            dst[Ipre, i, Ipost] = f(f(a1, a2), a3)
        end
        # process two edge points
        i = first(inds)
        a2 = dst[Ipre, i, Ipost]
        a3 = src[Ipre, i+1, Ipost]
        dst[Ipre, i, Ipost] = f(a2, a3)

        i = last(inds)
        a1 = src[Ipre, i-1, Ipost]
        a2 = dst[Ipre, i, Ipost]
        dst[Ipre, i, Ipost] = f(a1, a2)
    end
    return dst
end
