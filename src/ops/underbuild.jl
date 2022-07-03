"""
    underbuild(marker, mask; [dims])
    underbuild(marker, mask, se)

Reconstruction by dilation. This is an alias for [`mreconstruct`](@ref mreconstruct) with
`op=dilate`.

See also the in-place version [`underbuild!`](@ref), and the dual operator
[`overbuild`](@ref).
"""
underbuild(args...; kwargs...) = mreconstruct(dilate, args...; kwargs...)

"""
    underbuild!(out, marker, mask; [dims])
    underbuild!(out, marker, mask, se)

The in-place version of [`underbuild`](@ref) with output image `out` being modified in
place.
"""
underbuild!(args...; kwargs...) = mreconstruct!(dilate, args...; kwargs...)
