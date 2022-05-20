"""
    tophat(img; dims=coords_spatial(img), r=1)
    tophat(img, se)

Performs morphological top-hat transform for given image, i.e., `img - opening(img, se)`.

$(_docstring_se)

This _white_ top-hat transform can be used to extract small white elements and details from
an image. To extract black details, the _black_ top-hat transform, also known as bottom-hat
transform, [`bothat`](@ref) can be used.

# Examples

```jldoctest; setup = :(using ImageMorphology)
julia> img = falses(5, 5); img[1, 1] = true; img[3:5, 3:5] .= true; img
5×5 BitMatrix:
 1  0  0  0  0
 0  0  0  0  0
 0  0  1  1  1
 0  0  1  1  1
 0  0  1  1  1

julia> Int.(tophat(img))
5×5 Matrix{$Int}:
 1  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0

julia> Int.(tophat(img, strel_diamond(img))) # use diamond shape SE
5×5 Matrix{$Int}:
 1  0  0  0  0
 0  0  0  0  0
 0  0  1  0  0
 0  0  0  0  0
 0  0  0  0  0
```
"""
tophat(img; kwargs...) = tophat!(similar(img, maybe_floattype(eltype(img))), img, similar(img); kwargs...)
tophat(img, se) = tophat!(similar(img, maybe_floattype(eltype(img))), img, se, similar(img))

"""
    tophat!(out, img, buffer; dims=coords_spatial(img))
    tophat!(out, img, se, buffer)

The in-place version of [`tophat`](@ref) with input image `img` and output image `out`. The
intermediate erosion result is stored in `buffer`.
"""
function tophat!(out, img, buffer; dims=coords_spatial(img), r=nothing)
    return tophat!(out, img, strel_box(img, dims; r), buffer)
end
function tophat!(out, img, se, buffer)
    opening!(out, img, se, buffer)
    @. out = img - out
    return out
end
tophat!(out::AbstractArray{<:Color3}, img, se, buffer) = throw(ArgumentError("color image is not supported"))
