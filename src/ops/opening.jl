"""
    opening(img; dims=coords_spatial(img), r=1)
    opening(img, se)

Perform the morphological opening on `img`. The opening operation is defined as
[erosion](@ref erode) followed by a [dilation](@ref dilate): `dilate(erode(img, se), se)`.

$(_docstring_se)

# Examples

```jldoctest; setup = :(using ImageMorphology)
julia> img = trues(7,7); img[2, 2] = false; img[3:5, 3:5] .= false; img[4, 4] = true; img
7×7 BitMatrix:
 1  1  1  1  1  1  1
 1  0  1  1  1  1  1
 1  1  0  0  0  1  1
 1  1  0  1  0  1  1
 1  1  0  0  0  1  1
 1  1  1  1  1  1  1
 1  1  1  1  1  1  1

julia> opening(img)
7×7 BitMatrix:
 0  0  1  1  1  1  1
 0  0  1  1  1  1  1
 1  1  0  0  0  1  1
 1  1  0  0  0  1  1
 1  1  0  0  0  1  1
 1  1  1  1  1  1  1
 1  1  1  1  1  1  1

julia> opening(img, strel_diamond(img)) # use diamond shape SE
7×7 BitMatrix:
 1  1  1  1  1  1  1
 1  0  1  1  1  1  1
 1  1  0  0  0  1  1
 1  1  0  0  0  1  1
 1  1  0  0  0  1  1
 1  1  1  1  1  1  1
 1  1  1  1  1  1  1
```

## See also

- [`opening!`](@ref) is the in-place version of this function.
- [`closing`](@ref) is the dual operator of `opening` in the sense that
  `complement.(opening(img)) == closing(complement.(img))`.
"""
opening(img; kwargs...) = opening!(similar(img), img, similar(img); kwargs...)
opening(img, se) = opening!(similar(img), img, se, similar(img))

"""
    opening!(out, img, buffer; dims=coords_spatial(img))
    opening!(out, img, se, buffer)

The in-place version of [`opening`](@ref) with input image `img` and output image `out`. The
intermediate erosion result is stored in `buffer`.
"""
function opening!(out, img, buffer; dims=coords_spatial(img), r=nothing)
    return opening!(out, img, strel_box(img, dims; r), buffer)
end
function opening!(out, img, se, buffer)
    erode!(buffer, img, se)
    dilate!(out, buffer, se)
    return out
end
opening!(out::AbstractArray{<:Color3}, img, se, buffer) = throw(ArgumentError("color image is not supported"))
