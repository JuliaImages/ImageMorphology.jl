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
