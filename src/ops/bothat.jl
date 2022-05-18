"""
    bothat(img; dims=coords_spatial(img))
    bothat(img, se)

Performs morphological bottom-hat transform for given image, i.e., `closing(img, se) - img`.

`se` is the structuring element that defines the neighborhood of the image. See
[`strel`](@ref) for more details. If `se` is not specified, then it will use the
[`strel_box`](@ref) with an extra keyword `dims` to control the dimensions to filter.

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

julia> bothat(img)
7×7 BitMatrix:
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0
 0  0  0  1  0  0  1
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0

julia> bothat(img, strel_diamond(img)) # use diamond shape SE
7×7 BitMatrix:
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
bothat(img; dims=coords_spatial(img)) = bothat(img, strel_box(img, dims))
bothat(img, se) = bothat!(similar(img), img, se, similar(img))

"""
    bothat!(out, img, buffer; dims=coords_spatial(img))
    bothat!(out, img, se, buffer)

The in-place version of [`bothat`](@ref) with input image `img` and output image `out`. The
intermediate dilation result is stored in `buffer`.
"""
bothat!(out, img, buffer; dims=coords_spatial(img)) = bothat!(out, img, strel_box(img, dims), buffer)
function bothat!(out, img, se, buffer)
    closing!(out, img, se, buffer)
    @. out = out - img
    return out
end
