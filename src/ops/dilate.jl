"""
    dilate(img; dims=coords_spatial(img), r=1)
    dilate(img, se)

Perform a max-filter over the neighborhood of `img`, specified by structuring element `se`.

$(_docstring_se)

# Examples

```jldoctest; setup = :(using ImageMorphology)
julia> img = falses(5, 5); img[3, [2, 4]] .= true; img
5×5 BitMatrix:
 0  0  0  0  0
 0  0  0  0  0
 0  1  0  1  0
 0  0  0  0  0
 0  0  0  0  0

julia> dilate(img)
5×5 BitMatrix:
 0  0  0  0  0
 1  1  1  1  1
 1  1  1  1  1
 1  1  1  1  1
 0  0  0  0  0

julia> dilate(img; dims=1)
5×5 BitMatrix:
 0  0  0  0  0
 0  1  0  1  0
 0  1  0  1  0
 0  1  0  1  0
 0  0  0  0  0

julia> dilate(img, strel_diamond(img)) # use diamond shape SE
5×5 BitMatrix:
 0  0  0  0  0
 0  1  0  1  0
 1  1  1  1  1
 0  1  0  1  0
 0  0  0  0  0
```

## See also

- [`dilate!`](@ref) is the in-place version of this function
- [`erode`](@ref) is the dual operator of `dilate` in the sense that
  `complement.(dilate(img)) == erode(complement.(img))`.

!!! note "symmetricity"
    If `se` is symmetric with repsect to origin, i.e., `se[b] == se[-b]` for any `b`, then
    dilation becomes the Minkowski sum: A⊕B={a+b|a∈A, b∈B}.

"""
dilate(img::AbstractArray; kwargs...) = dilate!(similar(img), img; kwargs...)
dilate(img::AbstractArray, se::AbstractArray) = dilate!(similar(img), img, se)

"""
    dilate!(out, img; [dims])
    dilate!(out, img, se)

The in-place version of [`dilate`](@ref) with input image `img` and output image `out`.
"""
function dilate!(out, img; dims=coords_spatial(img), r=nothing)
    return dilate!(out, img, strel_box(img, dims; r))
end
dilate!(out, img, se::AbstractArray) = extreme_filter!(max, out, img, se)
dilate!(out::AbstractArray{<:Color3}, img, se::AbstractArray) = throw(ArgumentError("color image is not supported"))
