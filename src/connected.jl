"""
    label = label_components(A; [dims=coords_spatial(A)], [r=1], [bkg])
    label = label_components(A, se; [bkg])

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
[`component_lengths`](@ref), [`component_indices`](@ref), [`component_subscripts`](@ref),
[`component_centroids`](@ref) for basic properties of the labeled components.
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

function label_components!(out::AbstractArray{T}, A::AbstractArray, se; bkg=zero(eltype(A))) where {T<:Integer}
    axes(out) == axes(A) || throw_dmm(axes(out), axes(A))
    se = _maybe_build_symmetric_strel(se) # compat patch
    is_symmetric(se) || throw(ArgumentError("Non-symmetric structuring element is not supported yet"))
    upper_se, _ = strel_split(CartesianIndex, se)
    fill!(out, zero(T))
    sets = DisjointMinSets{T}()
    sizehint!(sets.parents, floor(Int, sqrt(length(A))))
    @inbounds for i in CartesianIndices(A)
        val = A[i]
        val == bkg && continue
        label = typemax(T)    # sentinel value
        for Δi in upper_se
            ii = i + Δi
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
function find_root_unsafe!(sets::DisjointMinSets, m::Int)
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

"`component_boxes(labeled_array)` -> an array of bounding boxes for each label, including the background label 0"
function component_boxes(img::AbstractArray{Int})
    nd = ndims(img)
    n = [Vector{Int}[fill(typemax(Int), nd), fill(typemin(Int), nd)] for i in 0:maximum(img)]
    s = CartesianIndices(size(img))
    for i in 1:length(img)
        vcur = s[i]
        vmin = n[img[i] + 1][1]
        vmax = n[img[i] + 1][2]
        for d in 1:nd
            vmin[d] = min(vmin[d], vcur[d])
            vmax[d] = max(vmax[d], vcur[d])
        end
    end
    return map(x -> map(y -> tuple(y...), x), n)
end

"`component_lengths(labeled_array)` -> an array of areas (2D), volumes (3D), etc. for each label, including the background label 0"
function component_lengths(img::AbstractArray{Int})
    n = zeros(Int, maximum(img) + 1)
    for i in 1:length(img)
        n[img[i] + 1] += 1
    end
    return n
end

"`component_indices(labeled_array)` -> an array of pixels for each label, including the background label 0"
function component_indices(img::AbstractArray{Int})
    n = [Int[] for i in 0:maximum(img)]
    for i in 1:length(img)
        push!(n[img[i] + 1], i)
    end
    return n
end

"`component_subscripts(labeled_array)` -> an array of pixels for each label, including the background label 0"
function component_subscripts(img::AbstractArray{Int})
    n = [Tuple[] for i in 0:maximum(img)]
    s = CartesianIndices(size(img))
    for i in 1:length(img)
        push!(n[img[i] + 1], s[i])
    end
    return n
end

"`component_centroids(labeled_array)` -> an array of centroids for each label, including the background label 0"
function component_centroids(img::AbstractArray{Int,N}) where {N}
    len = length(0:maximum(img))
    n = fill(zero(CartesianIndex{N}), len)
    counts = fill(0, len)
    @inbounds for I in CartesianIndices(size(img))
        v = img[I] + 1
        n[v] += I
        counts[v] += 1
    end
    return map(v -> n[v].I ./ counts[v], 1:len)
end
