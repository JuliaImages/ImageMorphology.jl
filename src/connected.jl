"""
    label = label_components(A; [dims=coords_spatial(A)], [r=1], [bkg], [periodic])
    label = label_components(A, se; [bkg], [periodic])

Find and label the connected components of array `A` where the connectivity is defined by
structuring element `se`. Each component is assigned a unique integer value as its label
with `0` representing the background defined by `bkg`.

$(_docstring_se)

# Examples

```jldoctest; setup=:(using ImageMorphology)
julia> A = [false true false true  false;
            true false false  true  true]
2×5 Matrix{Bool}:
 0  1  0  1  0
 1  0  0  1  1

julia> label_components(A) # default diamond shape C4 connectivity
2×5 Matrix{$Int}:
 0  2  0  3  0
 1  0  0  3  3

julia> label_components(A; dims=2) # only the rows are considered
2×5 Matrix{$Int}:
 0  2  0  3  0
 1  0  0  4  4

julia> label_components(A, strel_box((3, 3))) # box shape C8 connectivity
2×5 Matrix{$Int}:
 0  1  0  2  0
 1  0  0  2  2
```

The in-place version is [`label_components!`](@ref). See also [`component_boxes`](@ref),
[`component_lengths`](@ref), [`component_indices`](@ref), [`component_centroids`](@ref) for
basic properties of the labeled components.
"""
label_components(A; dims=coords_spatial(A), r=1, kwargs...) = label_components(A, strel_diamond(A, dims; r); kwargs...)
label_components(A, se; kwargs...) = label_components!(similar(A, Int), A, se; kwargs...)

"""
    label_components!(out, A; [dims], [r] [bkg])
    label_components!(out, A, se; [bkg])

The in-place version of [`label_components`](@ref).
"""
function label_components!(out, A; dims=coords_spatial(A), r=1, kwargs...)
    if !(axes(out) == axes(A))
        msg = "axes of input and output must match, got $(axes(out)) and $(axes(A))"
        hint = if eltype(A) <: Bool || eltype(A) <: CartesianIndex
            "The second argument seems to be a structuring element, it is expected to be the input array."
        else
            ""
        end
        throw(DimensionMismatch("$msg. $hint"))
    end
    return label_components!(out, A, strel_diamond(A, dims; r); kwargs...)
end

function label_components!(out::AbstractArray{T}, A::AbstractArray, se; bkg=zero(eltype(A)), periodic=false) where {T<:Integer}
    axes(out) == axes(A) || throw_dmm(axes(out), axes(A))
    se = _maybe_build_symmetric_strel(se) # compat patch
    is_symmetric(se) || throw(ArgumentError("Non-symmetric structuring element is not supported yet"))
    upper_se, _ = strel_split(CartesianIndex, se)
    fill!(out, zero(T))
    sets = DisjointMinSets{T}()
    sizehint!(sets.parents, floor(Int, sqrt(length(A))))
    # Define a function to compute indices for periodic boundary conditions
    periodic_index(i, Δi, sz) = mod1.(i .+ Δi .- 1, sz) .+ 1
    # Loop through all indices of A
    @inbounds for i in CartesianIndices(A)
        val = A[i]
        val == bkg && continue
        label = typemax(T)    # sentinel value
        for Δi in upper_se
            ii = i + Δi
            # If periodic, compute the periodic index
            if periodic
                ii = CartesianIndex(periodic_index(i, Δi, size(A)[k]) for k in 1:ndims(A))
            end
            checkbounds(Bool, A, ii) || continue
            if A[ii] == val
                newlabel = out[ii]
                label = if ((label == typemax(T)) | (label == newlabel))
                    newlabel
                else
                    union!(sets, label, newlabel)
                end
            end
        end
        if label == typemax(T)
            label = push!(sets)
        end
        out[i] = label
    end
    # Now parse sets to find the labels
    newlabel = minlabel(sets)
    @inbounds for i in eachindex(A, out)
        if A[i] != bkg
            out[i] = newlabel[find_root!(sets, out[i])]
        end
    end
    return out
end

function _maybe_build_symmetric_strel(iter)
    # NOTE(johnnychen94): label_components! (<v0.4) used to accept only the splited upper
    # half part as `se`. We should still support it with a deprecation warning.
    new_se = strel(CartesianIndex, strel(Bool, iter)) # remove the center position
    sym_se = union(.-new_se, new_se)
    if length(sym_se) == 2length(new_se)
        msg =
            "`label_components!` now requires a complete symmetric structuring element," *
            " but the provided structuring element is only the upper half." *
            " This is a deprecated behavior and will be removed in a future version."
        Base.depwarn(msg, :label_components!)
        return sym_se
    else
        return iter
    end
end

function throw_dmm(ax1, ax2)
    throw(DimensionMismatch("axes of input and output must match, got $ax1 and $ax2"))
end

# Copied directly from DataStructures.jl, but specialized
# to always make the parent be the smallest label
struct DisjointMinSets{T<:Integer}
    parents::Vector{T}

    DisjointMinSets{T}(n::Integer) where {T} = new([T(1):T(n);])
end
DisjointMinSets{T}() where {T<:Integer} = DisjointMinSets{T}(T(0))
DisjointMinSets() = DisjointMinSets{INt}()

Base.@propagate_inbounds function find_root!(sets::DisjointMinSets, m::Integer)
    p = sets.parents[m]   # don't use @inbounds here, it might not be safe
    @inbounds if sets.parents[p] != p
        sets.parents[m] = p = find_root_unsafe!(sets, p)
    end
    return p
end

# an unsafe variant of the above
function find_root_unsafe!(sets::DisjointMinSets, m::Integer)
    @inbounds p = sets.parents[m]
    @inbounds if sets.parents[p] != p
        sets.parents[m] = p = find_root_unsafe!(sets, p)
    end
    return p
end

Base.@propagate_inbounds function union!(sets::DisjointMinSets, m::Integer, n::Integer)
    mp = find_root!(sets, m)
    np = find_root!(sets, n)
    if mp < np
        sets.parents[np] = mp
        return mp
    elseif np < mp
        sets.parents[mp] = np
        return np
    end
    return mp
end

function Base.push!(sets::DisjointMinSets)
    m = length(sets.parents) + 1
    m >= typemax(eltype(sets.parents)) && error("labels exhausted, use a larger integer type")
    push!(sets.parents, m)
    return m
end

function minlabel(sets::DisjointMinSets)
    out = Vector{Int}(undef, length(sets.parents))
    k = 0
    @inbounds for i in 1:length(sets.parents)
        if sets.parents[i] == i
            k += 1
        end
        out[i] = k
    end
    return out
end

"""
    boxes = component_boxes(labeled_array)

Calculates the minimal bounding boxes for each label including the background label. The
labels can be computed by [`label_components`](@ref).

Each bounding box is represented as a `CartesianIndices`. `boxes` is shifted to 0-based
indexing vector so that background region is `boxes[0]`.

```jldoctest; setup=:(using ImageMorphology)
julia> A = [2 2 2 2 2; 1 1 1 0 1; 1 0 2 1 1; 1 1 2 2 2; 1 0 2 2 2]
5×5 Matrix{$Int}:
 2  2  2  2  2
 1  1  1  0  1
 1  0  2  1  1
 1  1  2  2  2
 1  0  2  2  2

julia> label = label_components(A) # four disjoint components
5×5 Matrix{$Int}:
 1  1  1  1  1
 2  2  2  0  4
 2  0  3  4  4
 2  2  3  3  3
 2  0  3  3  3

julia> boxes = component_boxes(label) # get bounding boxes of all regions
5-element OffsetArray(::Vector{CartesianIndices{2, Tuple{UnitRange{$Int}, UnitRange{$Int}}}}, 0:4) with eltype CartesianIndices{2, Tuple{UnitRange{$Int}, UnitRange{$Int}}} with indices 0:4:
 CartesianIndices((2:5, 2:4))
 CartesianIndices((1:1, 1:5))
 CartesianIndices((2:5, 1:3))
 CartesianIndices((3:5, 3:5))
 CartesianIndices((2:3, 4:5))

julia> A[boxes[1]] # crop the image region with label 1
1×5 Matrix{$Int}:
 2  2  2  2  2

julia> A[boxes[4]] # crop the image region with label 4
2×2 Matrix{$Int}:
 0  1
 1  1
```
"""
function component_boxes(A::AbstractArray{T,N}) where {T<:Integer,N}
    mn, mx = extrema(A)
    if !(mn == 0 || mn == 1)
        throw(ArgumentError("The input labeled array should contain background label `0` as the minimum value"))
    end
    boxes = if mn == 1
        OffsetArray(Matrix{CartesianIndex{N}}(undef, mx, 2), 0, 0)
    elseif mn == 0
        OffsetArray(Matrix{CartesianIndex{N}}(undef, 1 + mx, 2), -1, 0)
    end
    R = CartesianIndices(A)
    boxes[:, 1] .= Ref(last(R))
    boxes[:, 2] .= Ref(first(R))

    @inbounds for i in R
        label_idx = A[i]
        # this actually does a few unnecessary comparison, but it's relatively harmless as
        # the main computation is limited by memory
        boxes[label_idx, 1] = min(i, boxes[label_idx, 1])
        boxes[label_idx, 2] = max(i, boxes[label_idx, 2])
    end

    # convert to CartesianIndices
    return boxes = [boxes[i, 1]:boxes[i, 2] for i in axes(boxes, 1)]
end

"""
    counts = component_lengths(labeled_array)

Count the number of each labels in the input labeled array. `counts` is shifted to 0-based
indexing vector so that the number of background pixels is `counts[0]`.

```jldoctest; setup=:(using ImageMorphology)
julia> A = [2 2 2 2 2; 1 1 1 0 1; 1 0 2 1 1; 1 1 2 2 2; 1 0 2 2 2]
5×5 Matrix{$Int}:
 2  2  2  2  2
 1  1  1  0  1
 1  0  2  1  1
 1  1  2  2  2
 1  0  2  2  2

julia> label = label_components(A) # four disjoint components
5×5 Matrix{$Int}:
 1  1  1  1  1
 2  2  2  0  4
 2  0  3  4  4
 2  2  3  3  3
 2  0  3  3  3

julia> component_lengths(label)
5-element OffsetArray(::Vector{$Int}, 0:4) with eltype $Int with indices 0:4:
 3
 5
 7
 7
 3
```

For gray images, labels can be computed by [`label_components`](@ref).
"""
function component_lengths(A::AbstractArray{T}) where {T<:Integer}
    mn, mx = extrema(A)
    if !(mn == 0 || mn == 1)
        throw(ArgumentError("The input labeled array should contain background label `0` as the minimum value"))
    end
    counts = zeros(Int, 0:mx)
    for i in eachindex(A)
        counts[A[i]] += 1
    end
    return counts
end

"""
    indices = component_indices([T], labeled_array)

Get the indices of each label in the input labeled array. `indices` is shifted to 0-based
indexing vector so that the indices of background pixels is `indices[0]`.

The optional type `T` can be either `Int`/`IndexLinear()` or
`CartesianIndex`/`IndexCartesian()` that is used to specify the type of the indices. The
default choice is `IndexStyle(labeled_array)`.

```jldoctest; setup=:(using ImageMorphology)
julia> A = [2 2 2 2 2; 1 1 1 0 1; 1 0 2 1 1; 1 1 2 2 2; 1 0 2 2 2]
5×5 Matrix{$Int}:
 2  2  2  2  2
 1  1  1  0  1
 1  0  2  1  1
 1  1  2  2  2
 1  0  2  2  2

julia> label = label_components(A) # four disjoint components
5×5 Matrix{$Int}:
 1  1  1  1  1
 2  2  2  0  4
 2  0  3  4  4
 2  2  3  3  3
 2  0  3  3  3

julia> indices = component_indices(label)
5-element OffsetArray(::Vector{Vector{$Int}}, 0:4) with eltype Vector{$Int} with indices 0:4:
 [8, 10, 17]
 [1, 6, 11, 16, 21]
 [2, 3, 4, 5, 7, 9, 12]
 [13, 14, 15, 19, 20, 24, 25]
 [18, 22, 23]

julia> indices = component_indices(CartesianIndex, label)
5-element OffsetArray(::Vector{Vector{CartesianIndex{2}}}, 0:4) with eltype Vector{CartesianIndex{2}} with indices 0:4:
 [CartesianIndex(3, 2), CartesianIndex(5, 2), CartesianIndex(2, 4)]
 [CartesianIndex(1, 1), CartesianIndex(1, 2), CartesianIndex(1, 3), CartesianIndex(1, 4), CartesianIndex(1, 5)]
 [CartesianIndex(2, 1), CartesianIndex(3, 1), CartesianIndex(4, 1), CartesianIndex(5, 1), CartesianIndex(2, 2), CartesianIndex(4, 2), CartesianIndex(2, 3)]
 [CartesianIndex(3, 3), CartesianIndex(4, 3), CartesianIndex(5, 3), CartesianIndex(4, 4), CartesianIndex(5, 4), CartesianIndex(4, 5), CartesianIndex(5, 5)]
 [CartesianIndex(3, 4), CartesianIndex(2, 5), CartesianIndex(3, 5)]
```

For gray images, labels can be computed by [`label_components`](@ref).
"""
component_indices(A) = component_indices(IndexStyle(A), A)
component_indices(::Type{T}, A) where {T<:Integer} = component_indices(IndexLinear(), A)
component_indices(::Type{T}, A) where {T<:CartesianIndex} = component_indices(IndexCartesian(), A)
component_indices(::IndexCartesian, A) = _component_indices(A, CartesianIndices(A))
component_indices(::IndexLinear, A) = _component_indices(A, LinearIndices(A))

function _component_indices(A::AbstractArray{<:Integer,N}, R::AbstractArray{T,N}) where {T,N}
    mn, mx = extrema(A)
    if !(mn == 0 || mn == 1)
        throw(ArgumentError("The input labeled array should contain background label `0` as the minimum value"))
    end
    indices = OffsetVector([T[] for _ in 0:mx], -1)
    @inbounds for i in R
        push!(indices[A[i]], i)
    end
    return indices
end

"""
    centroids = component_centroids(labeled_array)

Compute the centroid of each label in the input labeled array. `centroids` is shifted to
0-based indexing vector so that the centroid of background pixels is `centroids[0]`.

The centroid of a finite set `X`, also known as geometric center, is calculated using
`sum(X)/length(X)`. For label `i`, all (Cartesian) indices of pixels with label `i` are used
to build the set `X`

```jldoctest; setup=:(using ImageMorphology)
julia> A = [2 2 2 2 2; 1 1 1 0 1; 1 0 2 1 1; 1 1 2 2 2; 1 0 2 2 2]
5×5 Matrix{$Int}:
 2  2  2  2  2
 1  1  1  0  1
 1  0  2  1  1
 1  1  2  2  2
 1  0  2  2  2

julia> label = label_components(A) # four disjoint components
5×5 Matrix{$Int}:
 1  1  1  1  1
 2  2  2  0  4
 2  0  3  4  4
 2  2  3  3  3
 2  0  3  3  3

julia> component_centroids(label)
5-element OffsetArray(::Vector{Tuple{Float64, Float64}}, 0:4) with eltype Tuple{Float64, Float64} with indices 0:4:
 (3.3333333333333335, 2.6666666666666665)
 (1.0, 3.0)
 (3.142857142857143, 1.5714285714285714)
 (4.285714285714286, 3.857142857142857)
 (2.6666666666666665, 4.666666666666667)
```

For gray images, labels can be computed by [`label_components`](@ref).
"""
function component_centroids(A::AbstractArray{<:Integer,N}) where {N}
    mn, mx = extrema(A)
    if !(mn == 0 || mn == 1)
        throw(ArgumentError("The input labeled array should contain background label `0` as the minimum value"))
    end
    n = fill(zero(CartesianIndex{N}), mn:mx)
    counts = fill(0, mn:mx)
    @inbounds for I in CartesianIndices(size(A))
        v = A[I]
        n[v] += I
        counts[v] += 1
    end
    return map((v, l) -> v.I ./ l, n, counts)
end

# deprecated
"`component_subscripts(labeled_array)` -> an array of pixels for each label, including the background label 0"
function component_subscripts(img::AbstractArray{Int})
    Base.depwarn("`component_subscripts` is deprecated. Use `component_indices(CartesianIndex, A)` instead.", :component_subscripts)
    n = [Tuple[] for i in 0:maximum(img)]
    s = CartesianIndices(size(img))
    for i in 1:length(img)
        push!(n[img[i] + 1], s[i])
    end
    return n
end

"""
    label_flatzones(img; [dims])
    label_flatzones(img; se)

Find the flat zones in image 'img'. Flat zones are defined as
connected (respect to se) voxels that have the same value.

The `dims` keyword is used to specify the dimension to process by constructing the box shape
structuring element [`strel_box(img; dims)`](@ref strel_box). For generic structuring
element, the half-size is expected to be either `0` or `1` along each dimension.

The output `labeled image` is an integer array.
"""
function label_flatzones(img; dims=coords_spatial(img))
    return label_flatzones(img, strel_box(img, dims))
end

function label_flatzones(img, se)
    return label_flatzones!(similar(img, Int), img, se)
end

function label_flatzones!(out, img; dims=coords_spatial(img))
    return label_flatzones!(out, img, strel_box(img, dims))
end

function label_flatzones!(out, img, se)
    return _generic_naive_labeling!(out, img, se)
end

"""
    label_lambdaflatzones(img; [dims])
    label_lambdaflatzones(img; se)

Find the lambda flat zones in image 'img'. Lambda Flat zones are defined as
connected (respect to se) voxels that have diff(value)<lambda.

The `dims` keyword is used to specify the dimension to process by constructing the box shape
structuring element [`strel_box(img; dims)`](@ref strel_box). For generic structuring
element, the half-size is expected to be either `0` or `1` along each dimension.

The output `labeled image` is an integer array.
"""

function label_lambdaflatzones(img, lambda; dims=coords_spatial(img))
    return label_lambdaflatzones(img, lambda, strel_box(img, dims))
end

function label_lambdaflatzones(img, lambda, se)
    return label_lambdaflatzones!(similar(img, Int), img, lambda, se)
end

function label_lambdaflatzones!(out, img, lambda; dims=coords_spatial(img))
    return label_lambdaflatzones!(out, img, lambda, strel_box(img, dims))
end

function label_lambdaflatzones!(out, img, lambda, se)
    return _generic_naive_labeling!(out, img, se, (x) -> true, (x, y) -> abs(x - y) <= lambda)
end

function _generic_naive_labeling!(out, img, se, can_be_labelled=(x) -> true, is_connected=(x, y) -> x == y)
    N = ndims(img)

    axes(out) == axes(img) || throw(DimensionMismatch("images should have the same axes"))

    se_size = strel_size(se)
    if length(se_size) != N
        msg = "the input structuring element is not for $N dimensional array, instead it is for $(length(se_size)) dimensional array"
        throw(DimensionMismatch(msg))
    end
    if !all(x -> in(x, (1, 3)), strel_size(se))
        msg = "structuring element with half-size larger than 1 is invalid"
        throw(DimensionMismatch(msg))
    end
    se = strel(CartesianIndex, se)

    fill!(out, 0)
    visited = similar(Array{Bool}, axes(img))
    fill!(visited, false)

    propagationfront = Queue{CartesianIndex{N}}() #Use queue to store propagation front

    R = CartesianIndices(axes(img))
    curr_label = 1
    for i in R
        # compute the "connected" zones of the image
        @inbounds if !visited[i] #candidate point aka not visited yet, start "connected" zones
            # mark point as visited
            @inbounds visited[i] = true
            @inbounds out[i] = curr_label
            enqueue!(propagationfront, i)
            while !isempty(propagationfront)
                #get indice from propagationfront
                idx_front = dequeue!(propagationfront)
                # get value at current point
                @inbounds centervalue = img[idx_front]
                if can_be_labelled(centervalue) # eg could be specialize for no background
                    #examine neighborhood points
                    for Δi in se # examine neighborhoods
                        ii = idx_front + Δi
                        if checkbounds(Bool, R, ii) #check that we are in the image
                            @inbounds nl_value = img[ii]
                            # is this pixel is connected ?
                            if is_connected(nl_value, centervalue)
                                # yes, propagate the "connected" zones to this pixel
                                @inbounds if !visited[ii]
                                    enqueue!(propagationfront, ii)
                                    # mark it as visited
                                    visited[ii] = true
                                    #propagate label
                                    out[ii] = curr_label
                                end
                            end
                        end
                    end #for all neighborValue
                end
            end # propagation front empty
            curr_label += 1 #increment label
            curr_label >= typemax(Int) && error("labels exhausted")
        end # test we can visit current point
    end # loop along idx image
    return out
end