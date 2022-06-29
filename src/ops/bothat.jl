"""
    bothat(img; dims=coords_spatial(img), r=1)
    bothat(img, se)

Performs morphological bottom-hat transform for given image, i.e., `closing(img, se) - img`.

$(_docstring_se)

This bottom-hat transform, also known as _black_ top-hat transform, can be used to extract
small black elements and details from an image. To extract white details, the _white_
top-hat transform [`tophat`](@ref) can be used.

# Examples

```jldoctest; setup = :(using ImageMorphology)
julia> img = falses(7, 7); img[3:5, 3:5] .= true; img[4, 6] = true; img[4, 4] = false; img
7×7 BitMatrix:
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0
 0  0  1  1  1  0  0
 0  0  1  0  1  1  0
 0  0  1  1  1  0  0
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0

julia> Int.(bothat(img))
7×7 Matrix{$Int}:
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0
 0  0  0  1  0  0  1
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0

julia> Int.(bothat(img, strel_diamond(img))) # use diamond shape SE
7×7 Matrix{$Int}:
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0
 0  0  0  1  0  0  0
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0
```

See also [`bothat!`](@ref) for the in-place version.
"""
bothat(img; kwargs...) = bothat!(similar(img, maybe_floattype(eltype(img))), img, similar(img); kwargs...)
bothat(img, se) = bothat!(similar(img, maybe_floattype(eltype(img))), img, se, similar(img))

"""
    bothat!(out, img, buffer; [dims], [r])
    bothat!(out, img, se, buffer)

The in-place version of [`bothat`](@ref) with input image `img` and output image `out`. The
intermediate dilation result is stored in `buffer`.
"""
function bothat!(out, img, buffer; dims=coords_spatial(img), r=nothing)
    return bothat!(out, img, strel_box(img, dims; r), buffer)
end
function bothat!(out, img, se, buffer)
    closing!(out, img, se, buffer)
    @. out = out - img
    return out
end
bothat!(out::AbstractArray{<:Color3}, img, se, buffer) = throw(ArgumentError("color image is not supported"))
