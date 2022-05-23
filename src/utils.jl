#=
    maybe_floattype(T)

To keep consistant with Base `diff` that Int array outputs Int array. It is sometimes
useful to only promote types for Bool and FixedPoint. In most of the time, `floattype`
should be the most reliable way.
=#
maybe_floattype(::Type{T}) where {T} = T
maybe_floattype(::Type{Bool}) = floattype(Bool)
maybe_floattype(::Type{T}) where {T<:FixedPoint} = floattype(T)
maybe_floattype(::Type{CT}) where {CT<:Color} = base_color_type(CT){maybe_floattype(eltype(CT))}

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

#=
A helper to eagerly check if particular function makes sense for `extreme_filter` semantics.
For instance, `max(::RGB, ::RGB)` is not well-defined and we should early throw errors so that
our users don't get encrypted error messages.
=#
function require_select_function(f, ::Type{T}) where {T}
    if !_is_select_function(f, T)
        hint = "does `f(x::T, y::T)` work as expected?"
        throw(ArgumentError("function `$f` is not a well-defined select function on type `$T`: $hint"))
    end
end
_is_select_function(f, ::Type{T}) where {T} = _is_select_function_trial(f, T)
function _is_select_function(f, ::Type{T}) where {T<:Real}
    f in (min, max) && return true
    return _is_select_function_trial(f, T)
end
function _is_select_function(f, ::Type{CT}) where {CT<:AbstractGray}
    f in (min, max) && return true
    return _is_select_function_trial(f, CT)
end
function _is_select_function(f, ::Type{CT}) where {CT<:Colorant}
    # min/max is not well-defined on generic color space
    f in (min, max) && return false
    return _is_select_function_trial(f, CT)
end
function _is_select_function_trial(f, ::Type{T}) where {T}
    # for generic case, just run a trial and see if it doesn't error
    v = zero(T)
    try
        f(v, v)
        return true
    catch
        return false
    end
end
