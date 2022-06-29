"""
    closing(img; dims=coords_spatial(img), r=1)
    closing(img, se)

Perform the morphological closing on `img`. The closing operation is defined as dilation
followed by an erosion: `erode(dilate(img, se), se)`.

$(_docstring_se)

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
closing(img; kwargs...) = closing!(similar(img), img, similar(img); kwargs...)
closing(img, se) = closing!(similar(img), img, se, similar(img))

"""
    closing!(out, img, buffer; [dims], [r])
    closing!(out, img, se, buffer)

The in-place version of [`closing`](@ref) with input image `img` and output image `out`. The
intermediate dilation result is stored in `buffer`.
"""
function closing!(out, img, buffer; dims=coords_spatial(img), r=nothing)
    return closing!(out, img, strel_box(img, dims; r), buffer)
end
function closing!(out, img, se, buffer)
    dilate!(buffer, img, se)
    erode!(out, buffer, se)
    return out
end
closing!(out::AbstractArray{<:Color3}, img, se, buffer) = throw(ArgumentError("color image is not supported"))
