"""
    overbuild(marker, mask; [dims])
    overbuild(marker, mask, se)

Reconstruction by erosion. This is an alias for [`mreconstruct`](@ref mreconstruct) with
`op=erode`.

See also the in-place version [`overbuild!`](@ref), and the dual operator
[`underbuild`](@ref).
"""
overbuild(args...; kwargs...) = mreconstruct(erode, args...; kwargs...)

"""
    overbuild!(out, marker, mask; [dims])
    overbuild!(out, marker, mask, se)

The in-place version of [`overbuild`](@ref) with output image `out` being modified in place.
"""
overbuild!(args...; kwargs...) = mreconstruct!(erode, args...; kwargs...)
