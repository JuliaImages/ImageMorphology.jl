"""
    morphogradient(img; dims=coords_spatial(img), r=1)
    morphogradient(img, se)

Calculate morphological (Beucher) gradient of the image, i.e., `dilate(img, se) - erode(img,
se)`.

$(_docstring_se)

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

julia> Int.(morphogradient(img))
7×7 Matrix{$Int}:
 0  0  0  0  0  0  0
 0  1  1  1  1  1  0
 0  1  1  1  1  1  0
 0  1  1  0  1  1  0
 0  1  1  1  1  1  0
 0  1  1  1  1  1  0
 0  0  0  0  0  0  0

julia> Int.(morphogradient(img, strel_diamond(img))) # use diamond shape SE
7×7 Matrix{$Int}:
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
function morphogradient(img; dims=coords_spatial(img), r=nothing)
    return morphogradient(img, strel_box(img, dims; r))
end
function morphogradient(img::AbstractArray{T}, se) where {T}
    require_symmetric_strel(se)
    buffer = similar(img)
    out = dilate!(similar(img, maybe_floattype(T)), img, se)
    buffer = erode!(similar(img, maybe_floattype(T)), img, se)
    return out .- buffer
end
morphogradient(img::AbstractArray{<:Color3}, se) = throw(ArgumentError("color image is not supported"))
