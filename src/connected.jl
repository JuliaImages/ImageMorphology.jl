import Base.push!  # for DisjointMinSets

"""
    label = label_components(A; bkg = zero(eltype(A)), dims=coords_spatial(A))
    label = label_components(A, connectivity; bkg = zero(eltype(A)))

Find the connected components in an array `A`. Components are defined as
connected voxels that all have the same value distinct from `bkg`, which
corresponds to the "background" component.

Specify connectivity in one of three ways:

- A list indicating which dimensions are used to determine
  connectivity. For example, `dims = (1,3)` would not test neighbors along
  dimension 2 for connectivity. This corresponds to just the nearest neighbors,
  i.e., default 4-connectivity in 2d and 6-connectivity in 3d.

- An iterable `connectivity` object with `CartesianIndex` elements encoding the
  displacement of each checked neighbor.

- A symmetric boolean array of the same dimensionality as `A`, of size 1 or 3
  along each dimension. Each entry in the array determines whether a given
  neighbor is used for connectivity analyses. For example, in two dimensions
  `connectivity = trues(3,3)` would include all pixels that touch the
  current one, even the corners.

The output `label` is an integer array, where `bkg` elements get a value of 0.

# Examples

```jldoctest; setup=:(using ImageMorphology)
julia> A = [true false false true  false;
            true false true  true  true]
2×5 Matrix{Bool}:
 1  0  0  1  0
 1  0  1  1  1

julia> label_components(A)
2×5 Matrix{$Int}:
 1  0  0  2  0
 1  0  2  2  2

julia> label_components(A; dims=2)
2×5 Matrix{$Int}:
 1  0  0  4  0
 2  0  3  3  3
```
With `dims=2`, entries in `A` are connected if they are in the same row, but
not if they are in the same column.
"""
label_components(A::AbstractArray; bkg=zero(eltype(A)), dims=coords_spatial(A)) =
    label_components(A, half_diamond(A, dims); bkg=bkg)
label_components(A::AbstractArray, connectivity::AbstractArray{Bool}; bkg=zero(eltype(A))) =
    label_components(A, half_pattern(A, connectivity); bkg=bkg)
label_components(A::AbstractArray, iter; bkg=zero(eltype(A))) =
    label_components!(similar(A, Int), A, iter; bkg=bkg) # Int is a safe choice

label_components!(out::AbstractArray{<:Integer}, A::AbstractArray; bkg=zero(eltype(A)), dims=coords_spatial(A)) =
    label_components!(out, A, half_diamond(A, dims); bkg=bkg)
label_components!(out::AbstractArray{<:Integer}, A::AbstractArray, connectivity::AbstractArray{Bool}; bkg=zero(eltype(A))) =
    label_components!(out, A, half_pattern(A, connectivity); bkg=bkg)
function label_components!(out::AbstractArray{T}, A::AbstractArray, iter; bkg=zero(eltype(A))) where T<:Integer
    axes(out) == axes(A) || throw_dmm(axes(out), axes(A))
    fill!(out, zero(T))
    sets = DisjointMinSets{T}()
    sizehint!(sets.parents, floor(Int, sqrt(length(A))))
    @inbounds for i in CartesianIndices(A)
        val = A[i]
        val == bkg && continue
        label = typemax(T)    # sentinel value
        for Δi in iter
            ii = i + Δi
            checkbounds(Bool, A, ii) || continue
            if A[ii] == val
                newlabel = out[ii]
                label = ((label == typemax(T)) | (label == newlabel)) ? newlabel : union!(sets, label, newlabel)
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

throw_dmm(ax1, ax2) = throw(DimensionMismatch("axes of input and output must match, got $ax1 and $ax2"))

function half_diamond(A::AbstractArray{T,N}, dims) where {T,N}
    offsets = CartesianIndex{N}[]
    for d = 1:N
        if d ∈ dims
            push!(offsets, CartesianIndex(ntuple(i -> i == d ? -1 : 0, N)))
        end
    end
    return (offsets...,)   # returning as a tuple allows specialization
    # return offsets
end
half_diamond(A::AbstractArray{T,N}, ::Colon) where {T,N} = half_diamond(A, 1:N)

function half_pattern(A::AbstractArray{T,N}, connectivity::AbstractArray{Bool}) where {T,N}
    all(in((1,3)), size(connectivity)) || throw(ArgumentError("connectivity must have size 1 or 3 in each dimension"))
    for d = 1:ndims(connectivity)
        size(connectivity, d) == 1 || reverse(connectivity; dims=d) == connectivity || throw(ArgumentError("connectivity must be symmetric"))
    end
    center = CartesianIndex(map(axes(connectivity)) do ax
        (first(ax) + last(ax)) ÷ 2
    end)
    offsets = CartesianIndex{N}[]
    for i in CartesianIndices(connectivity)
        i == center && break   # we only need the ones that come prior
        if connectivity[i]
            push!(offsets, i - center)
        end
    end
    # return offsets
    return (offsets...,)   # returning as a tuple allows specialization
end

# Copied directly from DataStructures.jl, but specialized
# to always make the parent be the smallest label
struct DisjointMinSets{T<:Integer}
    parents::Vector{T}

    DisjointMinSets{T}(n::Integer) where T = new([T(1):T(n);])
end
DisjointMinSets{T}() where T<:Integer = DisjointMinSets{T}(T(0))
DisjointMinSets() = DisjointMinSets{INt}()

Base.@propagate_inbounds function find_root!(sets::DisjointMinSets, m::Integer)
    p = sets.parents[m]   # don't use @inbounds here, it might not be safe
    @inbounds if sets.parents[p] != p
        sets.parents[m] = p = find_root_unsafe!(sets, p)
    end
    p
end

# an unsafe variant of the above
function find_root_unsafe!(sets::DisjointMinSets, m::Int)
    @inbounds p = sets.parents[m]
    @inbounds if sets.parents[p] != p
        sets.parents[m] = p = find_root_unsafe!(sets, p)
    end
    p
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
    mp
end

function push!(sets::DisjointMinSets)
    m = length(sets.parents) + 1
    m >= typemax(eltype(sets.parents)) && error("labels exhausted, use a larger integer type")
    push!(sets.parents, m)
    m
end

function minlabel(sets::DisjointMinSets)
    out = Vector{Int}(undef, length(sets.parents))
    k = 0
    @inbounds for i = 1:length(sets.parents)
        if sets.parents[i] == i
            k += 1
        end
        out[i] = k
    end
    out
end

"`component_boxes(labeled_array)` -> an array of bounding boxes for each label, including the background label 0"
function component_boxes(img::AbstractArray{Int})
    nd = ndims(img)
    n = [Vector{Int64}[ fill(typemax(Int64),nd), fill(typemin(Int64),nd) ]
            for i=0:maximum(img)]
    s = CartesianIndices(size(img))
    for i=1:length(img)
        vcur = s[i]
        vmin = n[img[i]+1][1]
        vmax = n[img[i]+1][2]
        for d=1:nd
            vmin[d] = min(vmin[d], vcur[d])
            vmax[d] = max(vmax[d], vcur[d])
        end
    end
    map(x->map(y->tuple(y...),x),n)
end

"`component_lengths(labeled_array)` -> an array of areas (2D), volumes (3D), etc. for each label, including the background label 0"
function component_lengths(img::AbstractArray{Int})
    n = zeros(Int64,maximum(img)+1)
    for i=1:length(img)
        n[img[i]+1]+=1
    end
    n
end

"`component_indices(labeled_array)` -> an array of pixels for each label, including the background label 0"
function component_indices(img::AbstractArray{Int})
    n = [Int64[] for i=0:maximum(img)]
    for i=1:length(img)
      push!(n[img[i]+1],i)
    end
    n
end

"`component_subscripts(labeled_array)` -> an array of pixels for each label, including the background label 0"
function component_subscripts(img::AbstractArray{Int})
    n = [Tuple[] for i=0:maximum(img)]
    s = CartesianIndices(size(img))
    for i=1:length(img)
      push!(n[img[i]+1],s[i])
    end
    n
end

"`component_centroids(labeled_array)` -> an array of centroids for each label, including the background label 0"
function component_centroids(img::AbstractArray{Int,N}) where N
    len = length(0:maximum(img))
    n = fill(zero(CartesianIndex{N}), len)
    counts = fill(0, len)
    @inbounds for I in CartesianIndices(size(img))
        v = img[I] + 1
        n[v] += I
        counts[v] += 1
    end
    map(v -> n[v].I ./ counts[v], 1:len)
end
