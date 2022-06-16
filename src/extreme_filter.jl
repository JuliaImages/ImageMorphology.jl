const MAX_OR_MIN = Union{typeof(max),typeof(min)}

"""
    extreme_filter(f, A; r=1, [dims]) -> out
    extreme_filter(f, A, Ω) -> out

Filter the array `A` using select function `f(x, y)` for each Ω-neighborhood. The name
"extreme" comes from the fact that typical select function `f` choice is `min` and `max`.

For each pixel `p` in `A`, the select function `f` is applied to its Ω-neighborhood
iteratively in a `f(...(f(f(A[p], A[p+Ω[1]]), A[p+Ω[2]]), ...)` manner. For instance, in the
1-dimensional case, `out[p] = f(f(A[p], A[p-1]), A[p+1])` for each `p` is the default
behavior.

The Ω-neighborhood is defined by the `dims` or `Ω` argument. The `r` and `dims` keywords
specifies the box shape neighborhood `Ω` using [`strel_box`](@ref). The `Ω` is also known as
structuring element (SE), it can be either displacement offsets or bool array mask, please
refer to [`strel`](@ref) for more details.

# Examples

```jldoctest extreme_filter; setup=:(using ImageMorphology)
julia> M = [4 6 5 3 4; 8 6 9 4 8; 7 8 4 9 6; 6 2 2 1 7; 1 6 5 2 6]
5×5 Matrix{$Int}:
 4  6  5  3  4
 8  6  9  4  8
 7  8  4  9  6
 6  2  2  1  7
 1  6  5  2  6

julia> extreme_filter(max, M) # max-filter using 4 direct neighbors along both dimensions
5×5 Matrix{$Int}:
 8  9  9  9  8
 8  9  9  9  9
 8  9  9  9  9
 8  8  9  9  9
 6  6  6  7  7

julia> extreme_filter(max, M; dims=1) # max-filter along the first dimension (column)
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

```jldoctest extreme_filter
julia> Ω_mask = centered(Bool[1 1 0; 1 1 0; 1 0 0]) # custom neighborhood in mask format
3×3 OffsetArray(::Matrix{Bool}, -1:1, -1:1) with eltype Bool with indices -1:1×-1:1:
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
extreme_filter(f, A; r=nothing, dims=coords_spatial(A)) = extreme_filter(f, A, strel_box(A, dims; r))
extreme_filter(f, A, Ω::AbstractArray) = extreme_filter!(f, similar(A), A, Ω)

"""
    extreme_filter!(f, out, A; [r], [dims])
    extreme_filter!(f, out, A, Ω)

The in-place version of [`extreme_filter`](@ref) where `out` is the output array that gets
modified.
"""
extreme_filter!(f, out, A; r=nothing, dims=coords_spatial(A)) = extreme_filter!(f, out, A, strel_box(A, dims; r))
function extreme_filter!(f, out, A, Ω)
    axes(out) == axes(A) || throw(DimensionMismatch("axes(out) must match axes(A)"))
    require_select_function(f, eltype(A))
    return _extreme_filter!(strel_type(Ω), f, out, A, Ω)
end

_extreme_filter!(::MorphologySE, f, out, A, Ω) = _extreme_filter_generic!(f, out, A, Ω)
_extreme_filter!(::SEDiamond, f, out, A, Ω) = _extreme_filter_diamond!(f, out, A, Ω)
function _extreme_filter!(::MorphologySE, f::MAX_OR_MIN, out, A::AbstractArray{T}, Ω) where {T<:Union{Gray{Bool},Bool}}
    # NOTE(johnnychen94): empirical choice based on benchmark results (intel i9-12900k)
    true_ratio = gray(sum(A) / length(A)) # this usually takes <2% of the time but gives a pretty good hint to do the decision
    true_ratio == 1 && (out .= true; return out)
    true_ratio == 0 && (out .= false; return out)
    use_bool = prod(strel_size(Ω)) > 9 || (f === max && true_ratio > 0.8) || (f === min && true_ratio < 0.2)
    if use_bool
        return _extreme_filter_bool!(f, out, A, Ω)
    else
        return _extreme_filter_generic!(f, out, A, Ω)
    end
end
function _extreme_filter!(::SEDiamond, f::MAX_OR_MIN, out, A::AbstractArray{T}, Ω) where {T<:Union{Gray{Bool},Bool}}
    rΩ = strel_size(Ω) .÷ 2
    if ndims(A) == 2 && all(rΩ[1] .== rΩ)
        # super fast implementation and wins in all cases
        return _extreme_filter_diamond!(f, out, A, Ω)
    end

    # NOTE(johnnychen94): empirical choice based on benchmark results (intel i9-12900k)
    true_ratio = gray(sum(A) / length(A)) # this usually takes <2% of the time but gives a pretty good hint to do the decision
    true_ratio == 1 && (out .= true; return out)
    true_ratio == 0 && (out .= false; return out)
    use_bool = (f === max && true_ratio > 0.4) || (f === min && true_ratio < 0.6)
    if use_bool
        return _extreme_filter_bool!(f, out, A, Ω)
    else
        return _extreme_filter_diamond!(f, out, A, Ω)
    end
end


###
# Implementation details
###

function _extreme_filter_generic!(f, out, A, Ω)
    @debug "call the generic `extreme_filter` implementation" fname = _extreme_filter_generic!
    Ω = strel(CartesianIndex, Ω)
    δ = CartesianIndex(strel_size(Ω) .÷ 2)

    R = CartesianIndices(A)
    R_inner = (first(R) + δ):(last(R) - δ)

    # NOTE(johnnychen94): delibrately duplicate the kernel codes here for inner loop and
    # boundary loop, because keeping the big function body allows Julia to do further
    # optimization and thus significantly improves the overall performance.
    @inbounds for p in R_inner
        # for interior points, boundary check is unnecessary
        s = A[p]
        for o in Ω
            v = A[p + o]
            s = f(s, v)
        end
        out[p] = s
    end
    @inbounds for p in EdgeIterator(R, R_inner)
        # for edge points, skip if the offset exceeds the boundary
        s = A[p]
        for o in Ω
            q = p + o
            checkbounds(Bool, A, q) || continue
            s = f(s, A[q])
        end
        out[p] = s
    end
    return out
end

# optimized implementation for SEDiamond -- a typical case of separable filter
function _extreme_filter_diamond!(f, out, A, Ω::SEDiamondArray{N}) where {N}
    rΩ = strel_size(Ω) .÷ 2
    if N == 2 && f isa MAX_OR_MIN && all(rΩ[1] .== rΩ)
        iter = rΩ[1]
        return _extreme_filter_diamond_2D!(f, out, A, iter)
    else
        return _extreme_filter_diamond_generic!(f, out, A, Ω)
    end
end

function _extreme_filter_diamond_generic!(f, out, A, Ω::SEDiamondArray)
    @debug "call the optimized `extreme_filter` implementation for SEDiamond SE" fname = _extreme_filter_diamond!
    rΩ = strel_size(Ω) .÷ 2
    r = maximum(rΩ)

    # To avoid the result affected by loop order, we need two arrays
    src = (out === A) || (r > 1) ? copy(A) : A
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
            Rpre = CartesianIndices(inds[1:(d - 1)])
            Rpost = CartesianIndices(inds[(d + 1):end])
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
        for i in (first(inds) + 1):(last(inds) - 1)
            a1 = src[Ipre, i - 1, Ipost]
            a2 = dst[Ipre, i, Ipost]
            a3 = src[Ipre, i + 1, Ipost]
            dst[Ipre, i, Ipost] = f(f(a1, a2), a3)
        end
        # process two edge points
        i = first(inds)
        a2 = dst[Ipre, i, Ipost]
        a3 = src[Ipre, i + 1, Ipost]
        dst[Ipre, i, Ipost] = f(a2, a3)

        i = last(inds)
        a1 = src[Ipre, i - 1, Ipost]
        a2 = dst[Ipre, i, Ipost]
        dst[Ipre, i, Ipost] = f(a1, a2)
    end
    return dst
end

function _extreme_filter_diamond_2D!(f::MAX_OR_MIN, out, A, iter)
    @debug "call the AVX-enabled `extreme_filter` implementation for 2D SEDiamond" fname = _extreme_filter_diamond_2D!
    actual_out = out

    # To avoid the result affected by loop order, we need two arrays
    src = (out === A) || (iter > 1) ? copy(A) : A

    # LoopVectorization currently doesn't understand Gray and N0f8 types, thus we
    # reinterpret to its raw data type UInt8
    src = rawview(channelview(src))
    out = rawview(channelview(out))

    Ω = strel(CartesianIndex, strel_diamond((3, 3)))
    δ = CartesianIndex(1, 1)
    R = CartesianIndices(src)
    R_inner = (first(R) + δ):(last(R) - δ)

    ox, oy = first(R).I

    # applying radius=r filter is equivalent to applying radius=1 filter iter times
    for i in 1:iter
        Rx, Ry = (ox + 1):(size(src, 1) + ox - 2), (oy + 1):(size(src, 2) + oy - 2)
        @turbo warn_check_args = false for y in Ry, x in Rx
            tmp = f(f(src[x, y], src[x - 1, y]), src[x + 1, y])
            out[x, y] = f(f(tmp, src[x, y - 1]), src[x, y + 1])
        end

        @inbounds for p in EdgeIterator(R, R_inner)
            # for edge points, skip if the offset exceeds the boundary
            s = src[p]
            for o in Ω
                q = p + o
                checkbounds(Bool, src, q) || continue
                s = f(s, src[q])
            end
            out[p] = s
        end

        if iter > 1 && i < iter
            src .= out
        end
    end
    return actual_out
end

# optimized implementation for Bool inputs with max/min select function
# 1) use &&, || instead of max, min
# 2) short-circuit the result to avoid unnecessary indexing and computation
function _extreme_filter_bool!(f, out, A::AbstractArray{Bool}, Ω)
    @debug "call the optimized max/min `extreme_filter` implementation for boolean array" fname = _extreme_filter_bool!
    Ω = strel(CartesianIndex, Ω)
    δ = CartesianIndex(strel_size(Ω) .÷ 2)

    R = CartesianIndices(A)
    R_inner = (first(R) + δ):(last(R) - δ)

    select = _fast_select(f)
    @inbounds for p in R_inner
        out[p] = select(A, p, Ω)
    end
    for p in EdgeIterator(R, R_inner)
        out[p] = select(A, p, Ω)
    end
    return out
end
function _extreme_filter_bool!(f, out, A::AbstractArray{Gray{Bool}}, Ω)
    return _extreme_filter_bool!(f, out, reinterpret(Bool, A), Ω)
end
_fast_select(::typeof(max)) = _maximum_fast
_fast_select(::typeof(min)) = _minimum_fast
Base.@propagate_inbounds function _maximum_fast(A::AbstractArray{Bool}, p, Ω)
    rst = A[p]
    for o in Ω
        # the complicated control flow ruins the SIMD and thus for small Ω, the performance
        # will be worse
        q = p + o
        @boundscheck checkbounds(Bool, A, q) || continue
        x = A[q]
        x && return true
        rst = rst || x
    end
    return rst
end
Base.@propagate_inbounds function _minimum_fast(A::AbstractArray{Bool}, p, Ω)
    rst = A[p]
    for o in Ω
        # the complicated control flow ruins the SIMD and thus for small Ω, the performance
        # will be worse
        q = p + o
        @boundscheck checkbounds(Bool, A, q) || continue
        x = A[q]
        x || return false
        rst = rst && x
    end
    return rst
end
