"""
    mgradient(img; dims=coords_spatial(img), r=1)
    mgradient(img, se)

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

julia> Int.(mgradient(img))
7×7 Matrix{$Int}:
 0  0  0  0  0  0  0
 0  1  1  1  1  1  0
 0  1  1  1  1  1  0
 0  1  1  0  1  1  0
 0  1  1  1  1  1  0
 0  1  1  1  1  1  0
 0  0  0  0  0  0  0

julia> Int.(mgradient(img, strel_diamond(img))) # use diamond shape SE
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

- [`mlaplace`](@ref) for the laplacian operator.
- `ImageBase.FiniteDiff` also provides a few finite difference operators, including `fdiff`,
  `fgradient`, etc.
"""
function mgradient(img; dims=coords_spatial(img), r=nothing)
    return mgradient(img, strel_box(img, dims; r))
end
function mgradient(img::AbstractArray{T}, se) where {T}
    require_symmetric_strel(se)
    buffer = similar(img)
    out = dilate!(similar(img, maybe_floattype(T)), img, se)
    buffer = erode!(similar(img, maybe_floattype(T)), img, se)
    return out .- buffer
end
mgradient(img::AbstractArray{<:Color3}, se) = throw(ArgumentError("color image is not supported"))
