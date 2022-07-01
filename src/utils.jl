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

#=
A helper to eagerly check if particular function makes sense for `extreme_filter` semantics.
For instance, `max(::RGB, ::RGB)` is not well-defined and we should early throw errors so that
our users don't get encrypted error messages.
=#
require_select_function(f, ::Type{T}) where {T} = require_select_function(f, T, T)
function require_select_function(f, ::Type{T1}, ::Type{T2}) where {T1,T2}
    if !_is_select_function(f, T1, T2)
        hint = "does `f(x::T1, y::T2)` work as expected?"
        throw(ArgumentError("function `$f` is not a well-defined select function on type `$T1` and `$T2`: $hint"))
    end
end
_is_select_function(f, ::Type{T}) where {T} = _is_select_function(f, T, T)

_is_select_function(f, ::Type{T1}, ::Type{T2}) where {T1,T2} = _is_select_function_trial(f, T1, T2)
function _is_select_function(f, ::Type{T1}, ::Type{T2}) where {T1<:Real,T2<:Real}
    f in (min, max) && return true
    return _is_select_function_trial(f, T1, T2)
end
function _is_select_function(f, ::Type{CT1}, ::Type{CT2}) where {CT1<:AbstractGray,CT2<:AbstractGray}
    f in (min, max) && return true
    return _is_select_function_trial(f, CT1, CT2)
end
function _is_select_function(f, ::Type{CT1}, ::Type{CT2}) where {CT1<:Colorant,CT2<:Colorant}
    # min/max is not well-defined on generic color space
    f in (min, max) && return false
    return _is_select_function_trial(f, CT1, CT2)
end
function _is_select_function_trial(f, ::Type{T1}, ::Type{T2}) where {T1,T2}
    # for generic case, just run a trial and see if it doesn't error
    try
        f(zero(T1), zero(T2))
        return true
    catch
        return false
    end
end
