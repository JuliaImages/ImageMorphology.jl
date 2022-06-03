"""
    upper, lower = strel_split([T], se)

Split a symmetric structuring element into its upper and lower half parts based on its
center point.

For each element `o` in `strel(CartesianIndex, upper)`, its negative `-o` is an element of
`strel(CartesianIndex, lower)`. This function is not the inverse of [`strel_chain`](@ref).

The splited non-symmetric SE parts will be represented as array of `T`, where `T` is either
a `Bool` or `CartesianIndex`. By default, `T = eltype(se)`.

```jldoctest strel_split; setup=:(using ImageMorphology; using ImageMorphology.StructuringElements)
julia> se = strel_diamond((3, 3))
3×3 SEDiamondArray{2, 2, UnitRange{$Int}, 0} with indices -1:1×-1:1:
 0  1  0
 1  1  1
 0  1  0

julia> upper, lower = strel_split(se);

julia> upper
3×3 OffsetArray(::Matrix{Bool}, -1:1, -1:1) with eltype Bool with indices -1:1×-1:1:
 0  1  0
 1  1  0
 0  0  0

julia> lower
3×3 OffsetArray(::Matrix{Bool}, -1:1, -1:1) with eltype Bool with indices -1:1×-1:1:
 0  0  0
 0  1  1
 0  1  0
```

If the `se` is represented as displacement offset array, then the splited result will also
be displacement offset array:

```jldoctest strel_split
julia> se = strel(CartesianIndex, se)
4-element Vector{CartesianIndex{2}}:
 CartesianIndex(0, -1)
 CartesianIndex(-1, 0)
 CartesianIndex(1, 0)
 CartesianIndex(0, 1)

julia> upper, lower = strel_split(se);

julia> upper
2-element Vector{CartesianIndex{2}}:
 CartesianIndex(0, -1)
 CartesianIndex(-1, 0)

julia> lower
2-element Vector{CartesianIndex{2}}:
 CartesianIndex(1, 0)
 CartesianIndex(0, 1)
```
"""
strel_split(se) = strel_split(eltype(se), se)
function strel_split(T, se)
    se = strel(Bool, se)
    require_symmetric_strel(se)
    R = LinearIndices(se)
    c = R[OffsetArrays.center(R)...]

    upper = copy(se)
    lower = copy(se)
    upper[(c + 1):end] .= false
    lower[begin:(c - 1)] .= false

    return strel(T, upper), strel(T, lower)
end

"""
    is_symmetric(se)

Check if a given structuring element array `se` is symmetric with respect to its center
pixel.

More formally, this checks if `mask[I] == mask[-I]` for any valid `I ∈
CartesianIndices(mask)` in the connectivity mask represetation `mask = strel(Bool, se)`.
"""
function is_symmetric(se::AbstractArray)
    # first check the axes, and then the values
    se = centered(strel(Bool, se))
    all(r -> first(r) == -last(r), axes(se)) || return false
    R = CartesianIndices(map(r -> 0:maximum(r), axes(se)))
    return all(R) do o
        @inbounds se[o] == se[-o]
    end
end
is_symmetric(se::SEBoxArray) = true
is_symmetric(se::SEDiamondArray) = true

#=
Some morphological operation only makes sense for symmetric structuring elements.
Here we provide a checker in spirit of Base.require_one_based_idnexing.
=#
require_symmetric_strel(se) = is_symmetric(se) || throw(ArgumentError("structuring element must be symmetric with respect to its center"))
