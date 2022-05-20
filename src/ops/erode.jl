"""
    out = erode(img; dims=coords_spatial(img), r=1)
    out = erode(img, se)

Perform a min-filter over the neighborhood of `img`, specified by structuring element `se`.

$(_docstring_se)

# Examples

```jldoctest; setup = :(using ImageMorphology)
julia> img = trues(5, 5); img[3, [2, 4]] .= false; img
5×5 BitMatrix:
 1  1  1  1  1
 1  1  1  1  1
 1  0  1  0  1
 1  1  1  1  1
 1  1  1  1  1

julia> erode(img)
5×5 BitMatrix:
 1  1  1  1  1
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 1  1  1  1  1

julia> erode(img; dims=1)
5×5 BitMatrix:
 1  1  1  1  1
 1  0  1  0  1
 1  0  1  0  1
 1  0  1  0  1
 1  1  1  1  1

julia> erode(img, strel_diamond(img)) # use diamond shape SE
5×5 BitMatrix:
 1  1  1  1  1
 1  0  1  0  1
 0  0  0  0  0
 1  0  1  0  1
 1  1  1  1  1
```

## See also

- [`erode!`](@ref) is the in-place version of this function
- [`dilate`](@ref) is the dual operator of `erode` in the sense that
  `complement.(dilate(img)) == erode(complement.(img))`.

!!! note "symmetricity"
    If `se` is symmetric with repsect to origin, i.e., `se[b] == se[-b]` for any `b`, then
    erosion becomes the Minkowski difference: A⊖B={a-b|a∈A, b∈B}.
"""
erode(img::AbstractArray; kwargs...) = erode!(similar(img), img; kwargs...)
erode(img::AbstractArray, se::AbstractArray) = erode!(similar(img), img, se)

"""
    erode!(out, img; [dims])
    erode!(out, img, se)

The in-place version of [`erode`](@ref) with input image `img` and output image `out`.
"""
function erode!(out, img; dims=coords_spatial(img), r=nothing)
    return erode!(out, img, strel_box(img, dims; r))
end
erode!(out, img, se::AbstractArray) = extreme_filter!(min, out, img, se)
erode!(out::AbstractArray{<:Color3}, img, se::AbstractArray) = throw(ArgumentError("color image is not supported"))
