"""
    closing(img; [dims])
    closing(img, se)

Perform the morphological closing on `img`. Mathematically, a closing operation is
[dilation](@ref dilate) followed by a [erosion](@ref erode):

```julia
closing(img, se) = erode(dilate(img, se), se)
```

`se` is the structuring element that defines the neighborhood of the image. See
[`strel`](@ref) for more details. If `se` is not specified, then it will use the
[`strel_box`](@ref) with an extra keyword `dims` to control the dimensions to filter.

# Examples

```jldoctest; setup = :(using ImageMorphology)
julia> img = falses(7,7); img[2, 2] = true; img[3:5, 3:5] .= true; img[4, 4] = false; img
7×7 BitMatrix:
 0  0  0  0  0  0  0
 0  1  0  0  0  0  0
 0  0  1  1  1  0  0
 0  0  1  0  1  0  0
 0  0  1  1  1  0  0
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0

julia> closing(img)
7×7 BitMatrix:
 1  1  0  0  0  0  0
 1  1  0  0  0  0  0
 0  0  1  1  1  0  0
 0  0  1  1  1  0  0
 0  0  1  1  1  0  0
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0

julia> closing(img, strel_diamond(img)) # # use diamond shape SE
7×7 BitMatrix:
 0  0  0  0  0  0  0
 0  1  0  0  0  0  0
 0  0  1  1  1  0  0
 0  0  1  1  1  0  0
 0  0  1  1  1  0  0
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0
```

## See also

- [`opening!`](@ref) is the in-place version of this function.
- [`closing`](@ref) is the dual operator of `opening` in the sense that
  `complement.(opening(img)) == closing(complement.(img))`.
"""
closing(img; dims=coords_spatial(img)) = closing!(similar(img), img, similar(img); dims)
closing(img, se) = closing!(similar(img), img, se, similar(img))

"""
    closing!(out, img, buffer; dims=coords_spatial(img))
    closing!(out, img, se, buffer)

The in-place version of [`closing`](@ref) with input image `img` and output image `out`. The
intermediate dilation result is stored in `buffer`.
"""
closing!(out, img, buffer; dims=coords_spatial(img)) = closing!(out, img, strel_box(img, dims), buffer)
function closing!(out, img, se, buffer)
    dilate!(buffer, img, se)
    erode!(out, buffer, se)
    return out
end
