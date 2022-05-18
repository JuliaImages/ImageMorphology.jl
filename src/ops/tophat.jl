"""
    tophat(img; dims=coords_spatial(img))
    tophat(img, se)

Performs morphological top-hat transform for given image, i.e., `img - opening(img, se)`.

`se` is the structuring element that defines the neighborhood of the image. See
[`strel`](@ref) for more details. If `se` is not specified, then it will use the
[`strel_box`](@ref) with an extra keyword `dims` to control the dimensions to filter.

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

julia> tophat(img)
5×5 BitMatrix:
 1  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0

julia> tophat(img, strel_diamond(img)) # use diamond shape SE
5×5 BitMatrix:
 1  0  0  0  0
 0  0  0  0  0
 0  0  1  0  0
 0  0  0  0  0
 0  0  0  0  0
```
"""
tophat(img; dims=coords_spatial(img)) = tophat(img, strel_box(img, dims))
tophat(img, se) = tophat!(similar(img), img, se, similar(img))

"""
    tophat!(out, img, buffer; dims=coords_spatial(img))
    tophat!(out, img, se, buffer)

The in-place version of [`tophat`](@ref) with input image `img` and output image `out`. The
intermediate erosion result is stored in `buffer`.
"""
tophat!(out, img, buffer; dims=coords_spatial(img)) = tophat!(out, img, strel_box(img, dims), buffer)
function tophat!(out, img, se, buffer)
    opening!(out, img, se, buffer)
    @. out = img - out
    return out
end
