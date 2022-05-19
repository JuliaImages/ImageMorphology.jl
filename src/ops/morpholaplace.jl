"""
    morpholaplace(img; dims=coords_spatial(img))
    morpholaplace(img, se)

Calculate morphological laplacian of the image.

The lapalacian operator is defined as `∇⁺A - ∇⁻A` where `∇⁺A` is the external gradient `A -
erode(A, se)` and `∇⁻A` is the internal gradient `dilate(A, se) - A`. Thus the laplacian is
`dilate(A, se) + erode(A, se) - 2A` in morphology sense.

`se` is the structuring element that defines the neighborhood of the image. See
[`strel`](@ref) for more details. If `se` is not specified, then it will use the
[`strel_box`](@ref) with an extra keyword `dims` to control the dimensions to filter.

# Examples

```jldoctest; setup = :(using ImageMorphology)
julia> img = falses(7, 7); img[3:5, 3:5] .= true; img[4, 4] = false; img
7×7 BitMatrix:
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0
 0  0  1  1  1  0  0
 0  0  1  0  1  0  0
 0  0  1  1  1  0  0
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0

julia> Int.(morpholaplace(img))
7×7 Matrix{$Int}:
 0  0   0   0   0  0  0
 0  1   1   1   1  1  0
 0  1  -1  -1  -1  1  0
 0  1  -1   1  -1  1  0
 0  1  -1  -1  -1  1  0
 0  1   1   1   1  1  0
 0  0   0   0   0  0  0

julia> Int.(morpholaplace(img, strel_diamond(img))) # use diamond shape SE
7×7 Matrix{$Int}:
 0  0   0   0   0  0  0
 0  0   1   1   1  0  0
 0  1  -1  -1  -1  1  0
 0  1  -1   1  -1  1  0
 0  1  -1  -1  -1  1  0
 0  0   1   1   1  0  0
 0  0   0   0   0  0  0
```

## See also

- [`morphogradient`](@ref) for the gradient operator.
- `ImageBase.FiniteDiff` also provides a few finite difference operators, including `fdiff`,
  `fgradient`, etc.
"""
morpholaplace(img; dims=coords_spatial(img)) = morpholaplace(img, strel_box(img, dims))
function morpholaplace(img::AbstractArray{T}, se) where {T}
    require_symmetric_strel(se)
    out = dilate!(similar(img, maybe_floattype(T)), img, se)
    buffer = erode!(similar(img, maybe_floattype(T)), img, se)
    @. out = out + buffer - 2img
    return out
end
