"""
    mlaplace(img; dims=coords_spatial(img), r=1)
    mlaplace(img, se)

Calculate morphological laplacian of the image.

The morphological lapalacian operator is defined as `∇⁺A - ∇⁻A` where `∇⁺A` is the external
gradient `A - erode(A, se)` and `∇⁻A` is the internal gradient `dilate(A, se) - A`. Thus is
`dilate(A, se) + erode(A, se) - 2A`.

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

julia> Int.(mlaplace(img))
7×7 Matrix{$Int}:
 0  0   0   0   0  0  0
 0  1   1   1   1  1  0
 0  1  -1  -1  -1  1  0
 0  1  -1   1  -1  1  0
 0  1  -1  -1  -1  1  0
 0  1   1   1   1  1  0
 0  0   0   0   0  0  0

julia> Int.(mlaplace(img, strel_diamond(img))) # use diamond shape SE
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

- [`mlaplace!`](@ref) is the in-place version of this function.
- [`mgradient`](@ref) for the gradient operator.
- `ImageBase.FiniteDiff` also provides a few finite difference operators, including `fdiff`,
  `fgradient`, etc.
"""
function mlaplace(img; dims=coords_spatial(img), r=nothing)
    return mlaplace(img, strel_box(img, dims; r))
end
function mlaplace(img::AbstractArray{T}, se) where {T}
    out = similar(img, maybe_floattype(T))
    buffer = similar(img, maybe_floattype(T))
    return mlaplace!(out, img, se, buffer)
end

"""
    mlaplace!(out, img, buffer; [dims], [r])
    mlaplace!(out, img, se, buffer)

The in-place version of [`mlaplace`](@ref) with input image `img` and output image `out`.
The intermediate erosion result is stored in `buffer`.
"""
function mlaplace!(out, img, buffer; dims=coords_spatial(img), r=nothing)
    return mlaplace!(out, img, strel_box(img, dims; r), buffer)
end
function mlaplace!(out, img, se, buffer)
    require_symmetric_strel(se)
    dilate!(out, img, se)
    erode!(buffer, img, se)
    @. out = out + buffer - 2img
    return out
end
mlaplace!(out::AbstractArray{<:Color3}, img, se, buffer=nothing) = throw(ArgumentError("color image is not supported"))
