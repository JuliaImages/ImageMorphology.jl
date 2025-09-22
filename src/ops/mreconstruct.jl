"""
    mreconstruct(op, marker, mask; [dims])
    mreconstruct(op, marker, mask, se)

Morphological reconstruction of `marker` image by operation `op`.

The `op` argument is either `erode` or `dilate`, indicating reconstruction by erosion or by
dilation. The `mask` argument has the same shape as `marker` and is used to restrict the
output value range.

The `dims` keyword is used to specify the dimension to process by constructing the box shape
structuring element [`strel_box(marker; dims)`](@ref strel_box). For generic structuring
element, the half-size is expected to be either `0` or `1` along each dimension.

By definition, the reconstruction is done by applying `marker = select.(op(marker; dims),
mask)` repeatly until reaching stability. For dilation `op, select = dilate, min` and for
erosion `op, select = erode, max`.

# Examples

```jldoctest; setup=:(using ImageMorphology)
julia> marker = [0 0 0 0 0; 0 9 0 0 0; 0 0 0 0 0; 0 0 0 5 0; 0 0 0 0 0; 0 9 0 0 0]
6×5 Matrix{$Int}:
 0  0  0  0  0
 0  9  0  0  0
 0  0  0  0  0
 0  0  0  5  0
 0  0  0  0  0
 0  9  0  0  0

julia> mask = [9 0 0 0 0; 0 8 7 1 0; 0 9 0 4 0; 0 0 0 4 0; 0 0 6 5 6; 0 0 9 8 9]
6×5 Matrix{$Int}:
 9  0  0  0  0
 0  8  7  1  0
 0  9  0  4  0
 0  0  0  4  0
 0  0  6  5  6
 0  0  9  8  9

julia> mreconstruct(dilate, marker, mask) # equivalent to underbuild(marker, mask)
6×5 Matrix{$Int}:
 8  0  0  0  0
 0  8  7  1  0
 0  8  0  4  0
 0  0  0  4  0
 0  0  4  4  4
 0  0  4  4  4
```

# See also

The inplace version of this function is [`mreconstruct!`](@ref). There are also aliases
[`underbuild`](@ref) for reconstruction by dilation and [`overbuild`](@ref) for
reconstruction by erosion.

# References

- [1] L. Vincent, “Morphological grayscale reconstruction in image analysis: applications
  and efficient algorithms,” IEEE Trans. on Image Process., vol. 2, no. 2, pp. 176–201, Apr.
  1993, doi: 10.1109/83.217222.
- [2] P. Soille, Morphological Image Analysis. Berlin, Heidelberg: Springer Berlin
  Heidelberg, 2004. doi: 10.1007/978-3-662-05088-0.
"""
function mreconstruct(op, marker, mask; dims=coords_spatial(mask))
    return mreconstruct(op, marker, mask, strel_box(marker, dims))
end
function mreconstruct(op, marker, mask, se)
    axes(marker) == axes(mask) || throw(DimensionMismatch("`marker` and `mask` should have the same axes"))
    # This ensures that we support mixed types, e.g. marker as Matrix{Gray{N0f8}} and mask as Matrix{Float64}
    OT = ImageCore.MosaicViews.promote_wrapped_type(eltype(marker), eltype(mask))
    return mreconstruct!(op, similar(mask, OT), marker, mask, se)
end

"""
    mreconstruct!(op, out, marker, mask; [dims])

The in-place version of morphological reconstruction [`mreconstruct`](@ref).
"""
function mreconstruct!(op, out, marker, mask; dims=coords_spatial(mask))
    return mreconstruct!(op, out, marker, mask, strel_box(marker, dims))
end
function mreconstruct!(op, out, marker, mask, se)
    return _mreconstruct!(_prepare_reconstruct_ops(op), out, marker, mask, se)
end

_prepare_reconstruct_ops(::op) where {op} = error("operation `$(op.instance)` is not supported for `mreconstruct`")
_prepare_reconstruct_ops(::typeof(dilate)) = (max, min)
_prepare_reconstruct_ops(::typeof(erode)) = (min, max)
@inline _should_push(::typeof(max)) = <
@inline _should_push(::typeof(min)) = >

# Single-CPU version, references are [1] and section 6.2 of [2]
function _mreconstruct!((select_se, select_marker), out, marker, mask, se)
    N = ndims(marker)

    axes(out) == axes(marker) == axes(mask) || throw(DimensionMismatch("images should have the same axes"))
    require_select_function(select_se, eltype(mask), eltype(marker))
    require_select_function(select_marker, eltype(mask), eltype(marker))

    # For generic structuring element, the half-size should to be either `0` or `1` along
    # each dimension. See also section 6.2.3 of [2].
    se_size = strel_size(se)
    if length(se_size) != N
        msg = "the input structuring element is not for $N dimensional array, instead it is for $(length(se_size)) dimensional array"
        throw(DimensionMismatch(msg))
    end
    if !all(x -> in(x, (1, 3)), strel_size(se))
        # center cropping the input se using [-1:1, -1:1, ...]
        inds_str = join(ntuple(_ -> 3, N), "×")
        @warn "structuring element with half-size larger than 1 is invalid, only the center $inds_str values are used"
        inds = ntuple(i -> (-1:1), N)
        se = centered(strel(Bool, strel(CartesianIndex, se))[inds...])
    end
    se = strel(CartesianIndex, se)

    # The original algorithm assume that the marker < mask for erosion and marker > mask for
    # dilation. This assertion is not always true in practice.
    @. out = select_marker(marker, mask)

    queue = Queue{CartesianIndex{N}}()
    se = strel(CartesianIndex, se)
    upper_se, lower_se = strel_split(se)

    # NOTE we could pad the array with -Inf to speedup forward/backward scan
    # forward scan
    R = CartesianIndices(axes(marker))
    for i in R
        @inbounds curr_val = out[i]
        for Δi in upper_se # examine neighborhoods
            ii = i + Δi
            if checkbounds(Bool, R, ii) #check that we are in the image
                @inbounds curr_val = select_se(curr_val, out[ii])
            end
        end
        @inbounds out[i] = select_marker(curr_val, mask[i])
    end
    # backward scan
    should_push = _should_push(select_se)
    for i in reverse(R)
        @inbounds curr_val = out[i]
        for Δi in lower_se # examine neighborhoods
            ii = i + Δi
            if checkbounds(Bool, R, ii) #check that we are in the image
                @inbounds curr_val = select_se(curr_val, out[ii])
            end
        end
        @inbounds out[i] = select_marker(curr_val, mask[i])
        for Δi in lower_se # examine neighborhoods
            ii = i + Δi
            if checkbounds(Bool, R, ii) #check that we are in the image
                @inbounds if should_push(out[ii], out[i]) && should_push(out[ii], mask[ii])
                    push!(queue, i)
                end
            end
        end
    end
    # Loop until all pixel have been examined
    while !isempty(queue)
        curr_idx = popfirst!(queue)
        for Δi in se # examine neighborhoods
            ii = curr_idx + Δi
            if checkbounds(Bool, R, ii) #check that we are in the image
                @inbounds if should_push(out[ii], out[curr_idx]) && mask[ii] != out[ii]
                    out[ii] = select_marker(out[curr_idx], mask[ii])
                    push!(queue, ii)
                end
            end
        end
    end
    return out
end
