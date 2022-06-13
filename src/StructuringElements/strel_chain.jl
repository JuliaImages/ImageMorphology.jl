struct SEChain{N} <: MorphologySE{N}
    data::Vector{<:AbstractArray{Bool}}
end

struct SEChainArray{N,AT<:AbstractArray{Bool,N}} <: MorphologySEArray{N}
    data::Vector{<:AbstractArray{Bool}}
    _data::AT # pre-calculated chained data as mask
    function SEChainArray{N}(data::Vector{<:AbstractArray{Bool}}) where {N}
        # TODO(johnnychen94): defer the calculation of chained data until it is needed since
        # this is mainly used for visualization purposes
        _data = _strel_chain(data)
        return new{N,typeof(_data)}(data, _data)
    end
end
function SEChainArray(data::Vector{<:AbstractArray{Bool}})
    data = convert(Vector, data)
    N = maximum(ndims.(data))
    return SEChainArray{N}(data)
end
SEChainArray(data::Tuple) = SEChainArray(collect(data))
strel_type(A::SEChainArray{N}) where {N} = SEChain{N}(A.data)
strel_size(A::SEChainArray) = size(A._data)

Base.axes(A::SEChainArray) = axes(A._data)
Base.size(A::SEChainArray) = size(A._data)
Base.@propagate_inbounds Base.getindex(A::SEChainArray, inds::Int...) = getindex(A._data, inds...)

"""
    strel_chain(A, B, ...)
    strel_chain(As)

For structuring elements of the same dimensions, chain them together to build a bigger one.

The output dimension is the same as the inputs dimensions. See also [strel_product](@ref
strel_product) that cartesian producting each SE.

!!! note "structuring element decomposition"
    For some morphological operations `f` such as dilation and erosion, if `se` can be
    decomposed into smaller ones `se1, se2, ..., seN`, then `f(img, se)` is equivalent to
    `f(...f(f(img, se1), se2), ..., seN)`. Because applying `f` to smaller SEs is more
    efficient than to the original big one, this trick is used widely in image morphology.

```jldoctest; setup=:(using ImageMorphology; using ImageMorphology.StructuringElements)
julia> img = rand(512, 512);

julia> se1, se2 = [centered(rand(Bool, 3, 3)) for _ in 1:2];

julia> se = strel_chain(se1, se2);

julia> out_se = dilate(img, se);

julia> out_pipe = dilate(dilate(img, se1), se2);

julia> out_se[2:end-1, 2:end-1] == out_pipe[2:end-1, 2:end-1] # except for the boundary
true
```
"""
strel_chain(se, se_rest...) = strel_chain([se, se_rest...])
strel_chain(se_list::Vector{<:AbstractArray{T,N}}) where {T,N} = SEChainArray{N}(se_list)
strel_chain(se_list::Tuple) = SEChainArray(se_list)
strel_chain(se) = se

function _strel_chain(data::AbstractVector{<:AbstractArray})
    isempty(data) && throw(ArgumentError("data cannot be empty"))
    data = strel.(CartesianIndex, data)
    out = first(data)
    for i in axes(data, 1)[2:end]
        # TODO: preallocating the output can reduce a few more
        out = _simple_dilate(out, data[i])
    end
    out = strel(Bool, strel(CartesianIndex, out))
    return out
end

# a simple dilate function that automatically extends the boundary and dimension
function _simple_dilate(A::AbstractArray{T,N}, B::AbstractArray{T,N}) where {T,N}
    r = strel_size(B) .÷ 2
    sz = max.(strel_size(A), strel_size(B))
    out_sz = @. sz + 2r
    out = centered(falses(out_sz))

    R = strel(CartesianIndex, A)
    offsets = strel(CartesianIndex, B)
    i = zero(eltype(R))
    out[i] = true
    for o in offsets
        out[i + o] = true
    end
    for i in R
        out[i] = true
        for o in offsets
            out[i + o] = true
        end
    end

    # remove unnecessary zero boundaries
    return strel(CartesianIndex, out)
end

"""
    strel_product(A, B, ...)
    strel_product(se_list)

Cartesian product of multiple structuring elements; the output dimension `ndims(out) ==
sum(ndims, se_list)`.

See also [`strel_chain`](@ref) that chains SEs in the same dimension.

```jldoctest; setup=:(using ImageMorphology; using ImageMorphology.StructuringElements)
julia> strel_product(strel_diamond((5, 5)), centered(Bool[1, 1, 1]))
5×5×3 SEChainArray{3, OffsetArrays.OffsetArray{Bool, 3, BitArray{3}}} with indices -2:2×-2:2×-1:1:
[:, :, -1] =
 0  0  1  0  0
 0  1  1  1  0
 1  1  1  1  1
 0  1  1  1  0
 0  0  1  0  0

[:, :, 0] =
 0  0  1  0  0
 0  1  1  1  0
 1  1  1  1  1
 0  1  1  1  0
 0  0  1  0  0

[:, :, 1] =
 0  0  1  0  0
 0  1  1  1  0
 1  1  1  1  1
 0  1  1  1  0
 0  0  1  0  0
```
"""
strel_product(se, se_rest...) = strel_product([se, se_rest...])

# TODO(johnnychen94): fix the type instability if it really matters in practice
function strel_product(se_list)
    N = sum(ndims, se_list)
    new_se_list = Array{Bool,N}[]
    ni = 0
    for se in se_list
        pre_ones = ntuple(_ -> one(Int), max(0, ni))
        post_ones = ntuple(_ -> one(Int), max(0, N - ndims(se) - ni))
        new_se = reshape(se, pre_ones..., size(se)..., post_ones...)
        push!(new_se_list, convert(Array{Bool,N}, new_se))
        ni += ndims(se)
    end
    return SEChainArray{N}(map(centered, new_se_list))
end
