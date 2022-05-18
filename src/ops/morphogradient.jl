"""
    morphogradient(img; dims=coords_spatial(img))
    morphogradient(img, se)

Calculate morphological (Beucher) gradient of the image, i.e., `dilate(img, se) - erode(img,
se)`.

`se` is the structuring element that defines the neighborhood of the image. See
[`strel`](@ref) for more details. If `se` is not specified, then it will use the
[`strel_box`](@ref) with an extra keyword `dims` to control the dimensions to filter.

# Examples

```jldoctest; setup = :(using ImageMorphology)
julia> img = falses(7, 7); img[3:5, 3:5] .= true; img
7×7 BitMatrix:
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0
 0  0  1  1  1  0  0
 0  0  1  1  1  0  0
 0  0  1  1  1  0  0
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0

julia> BitArray(morphogradient(img))
7×7 BitMatrix:
 0  0  0  0  0  0  0
 0  1  1  1  1  1  0
 0  1  1  1  1  1  0
 0  1  1  0  1  1  0
 0  1  1  1  1  1  0
 0  1  1  1  1  1  0
 0  0  0  0  0  0  0

julia> BitArray(morphogradient(img, strel_diamond(img))) # use diamond shape SE
7×7 BitMatrix:
 0  0  0  0  0  0  0
 0  0  1  1  1  0  0
 0  1  1  1  1  1  0
 0  1  1  0  1  1  0
 0  1  1  1  1  1  0
 0  0  1  1  1  0  0
 0  0  0  0  0  0  0
```

## See also

- [`morpholaplace`](@ref) for the laplacian operator.
- `ImageBase.FiniteDiff` also provides a few finite difference operators, including `fdiff`,
  `fgradient`, etc.
"""
morphogradient(img; dims=coords_spatial(img)) = morphogradient(img, strel_box(img, dims))
function morphogradient(img::AbstractArray{T}, se) where {T}
    buffer = similar(img)
    out = dilate!(similar(img, maybe_floattype(T)), img, se)
    buffer = erode!(similar(img, maybe_floattype(T)), img, se)
    return out .- buffer
end
