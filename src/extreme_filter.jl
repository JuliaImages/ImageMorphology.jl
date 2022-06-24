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

# Shift vector of size N up or down by 1 accorded to dir, pad with v
const SafeArrayTypes{T,N,NA} = Union{Array{T,N},Base.SubArray{T,N,Array{T,NA}}}
function _unsafe_padded_copyto!(dest::SafeArrayTypes{T,1}, src::SafeArrayTypes{T,1}, dir, N, v) where {T}
    # ptr level optimized implementation for Real types
    # short-circuit all check so unsafe
    if dir
        unsafe_store!(pointer(dest), v)
        unsafe_copyto!(pointer(dest, 2), pointer(src, 1), N - 1)
    else
        unsafe_copyto!(pointer(dest, 1), pointer(src, 2), N - 1)
        unsafe_store!(pointer(dest), v, N)
    end
    return dest
end
function _unsafe_padded_copyto!(dest::AbstractVector, src::AbstractVector, dir, N, v)
    if dir
        dest[begin] = v
        copyto!(dest, 2, src, 1, N - 1)
    else
        dest[end] = v
        copyto!(dest, 1, src, 2, N - 1)
    end
    return dest
end

# ptr level optimized implementation for Real types
# short-circuit all check
# call LoopVectorization directly
function _unsafe_shift_arith!(
    f::MAX_OR_MIN,
    out::AbstractVector,
    tmp::AbstractVector,
    tmp2::AbstractVector,
    A::AbstractVector,
) #tmp external to reuse external allocation
    if f === min
        padd = typemax(eltype(out))
    else
        padd = typemin(eltype(out))
    end
    N = length(out)
    #in  src  = [1,2,3,4], padd, N=4
    #out tmp  = [padd,1,2,3]
    #out tmp2 = [1,2,3,padd]
    _unsafe_padded_copyto!(tmp, A, true, N, padd)
    _unsafe_padded_copyto!(tmp2, A, false, N, padd)
    @turbo warn_check_args = false for i in eachindex(out, A, tmp, tmp2)
        out[i] = f(f(A[i], tmp[i]), tmp2[i])
    end
    return out
end

# optimized `extreme_filter` implementation for SEDiamond SE and 2D images
# Similar approach could be found in
# Žlaus, Danijel & Mongus, Domen. (2018). In-place SIMD Accelerated Mathematical Morphology. 76-79. 10.1145/3206157.3206176.
# see https://www.researchgate.net/publication/325480366_In-place_SIMD_Accelerated_Mathematical_Morphology
function _extreme_filter_diamond_2D!(f::MAX_OR_MIN, out, A, iter)
    @debug "call the AVX-enabled `extreme_filter` implementation for 2D SEDiamond" fname = _extreme_filter_diamond_2D!
    out_actual = out

    out = OffsetArrays.no_offset_view(out)
    A = OffsetArrays.no_offset_view(A)

    # To avoid the result affected by loop order, we need two arrays
    src = (out === A) || (iter > 1) ? copy(A) : A
    # NOTE(johnnychen94): we don't need to do `out .= src` here because it is write-only;
    # all values are generated from the read-only `src`.

    # LoopVectorization currently doesn't understand Gray and N0f8 types, thus we
    # reinterpret to its raw data type UInt8
    src = rawview(channelview(src))
    out = rawview(channelview(out))

    # creating temporaries
    ySize, xSize = size(A)
    tmp = similar(out, eltype(out), ySize)
    tmp2 = similar(tmp)
    tmp3 = similar(tmp)

    # applying radius=r filter is equivalent to applying radius=1 filter r times
    for i in 1:iter
        # NOTE(johnnychen94): this implementation is essentially equivalent to the
        # following loop. Here we buffer the getindex to further improve the performance
        # by explicitly introducing SIMD buffers.

        # for y in 2:(size(src, 2) - 1), x in 2:(size(src, 1) - 1)
        #     tmp = f(f(src[x, y], src[x - 1, y]), src[x + 1, y])
        #     out[x, y] = f(f(tmp, src[x, y - 1]), src[x, y + 1])
        # end

        #compute first edge column
        #translate to clipped connection, x neighborhood of interest, ? neighborhood we don't care, . center
        # x ?
        # . x
        # x ?
        viewprevious = view(src, :, 1)
        # dilate/erode col 1
        _unsafe_shift_arith!(f, tmp2, tmp, tmp3, viewprevious)
        viewnext = view(src, :, 2)
        viewout = view(out, :, 1)
        # inf/sup between dilate/erode col 1 0 and col 2
        LoopVectorization.vmap!(f, viewout, viewnext, tmp2)
        #next->current
        viewcurrent = view(src, :, 2)
        for c in 3:xSize
            # viewprevious c-2
            # viewcurrent  c-1
            # viewnext     c
            # ? x ?
            # x . x
            # ? x ?
            viewout = view(out, :, c - 1)
            viewnext = view(src, :, c)
            # dilate(x-1)/erode(x-1)
            _unsafe_shift_arith!(f, tmp2, tmp, tmp3, viewcurrent)
            @turbo warn_check_args = false for i in eachindex(viewout, viewprevious, tmp2)
                #sup(x-2,dilate(x-1)),inf(x-2,erode(x-1))
                #sup(sup(x-2,dilate(x-1),x) || inf(inf(x-2,erode(x-1),x)
                viewout[i] = f(f(viewprevious[i], tmp2[i]), viewnext[i])
            end
            #current->previous
            viewprevious = view(src, :, c - 1)
            #next->current
            viewcurrent = view(src, :, c)
        end
        #end last column
        #translate to clipped connection
        # ? x
        # x .
        # ? x
        # dilate/erode col x
        viewout = view(out, :, xSize)
        _unsafe_shift_arith!(f, tmp2, tmp, tmp3, viewcurrent)
        LoopVectorization.vmap!(f, viewout, tmp2, viewprevious)
        if iter > 1 && i < iter
            src .= out
        end
    end

    return out_actual
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
