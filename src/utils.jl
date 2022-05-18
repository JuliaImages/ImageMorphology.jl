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
