"""
    mgradient(img; mode=:beucher, dims=coords_spatial(img), r=1)
    mgradient(img, se; mode=:beucher)

Calculate morphological gradient of the image using given mode.

There are three widely used modes[1]:

- `:beucher`: the default mode. It calculates the arithmetic difference between the dilation
  and the erosion -- `dilate(img, se) - erode(img, se)`.
- `:internal`: also known as _half-gradient by erosion_. It calculates the arithmetic
  difference between the original image and its erosion -- `img - erode(img, se)`.
- `:external`: also known as _half-gradient by dilation_. It calculates the arithmetic
  difference between dilation and the original image -- `dilate(img, se) - se`.

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

julia> Int.(mgradient(img)) # default mode :beucher always creates a two-pixel wide boundary
7×7 Matrix{$Int}:
 0  0  0  0  0  0  0
 0  1  1  1  1  1  0
 0  1  1  1  1  1  0
 0  1  1  0  1  1  0
 0  1  1  1  1  1  0
 0  1  1  1  1  1  0
 0  0  0  0  0  0  0

julia> Int.(mgradient(img; mode=:internal)) # half-gradient -- the boundary is internal to original image
7×7 Matrix{$Int}:
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0
 0  0  1  1  1  0  0
 0  0  1  0  1  0  0
 0  0  1  1  1  0  0
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0

julia> Int.(mgradient(img; mode=:external)) # half-gradient -- the boundary is external to original image
7×7 Matrix{$Int}:
 0  0  0  0  0  0  0
 0  1  1  1  1  1  0
 0  1  0  0  0  1  0
 0  1  0  0  0  1  0
 0  1  0  0  0  1  0
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

The beucher operator is a self-complementary operator in the sense that `mgradient(img, se;
mode=:beucher) == mgradient(complement.(img), se; mode=:beucher)`. When `r>1`, it is usually
called _thick gradient_. If a line segment is used as `se`, then the gradient becomes the
_directional gradient_.

## See also

- [`mgradient!`](@ref) is the in-place version of this function.
- [`mlaplace`](@ref) for the laplacian operator.
- `ImageBase.FiniteDiff` also provides a few finite difference operators, including `fdiff`,
  `fgradient`, etc.

## References

- [1] Rivest, Jean-Francois, Pierre Soille, and Serge Beucher. "Morphological gradients."
  Journal of Electronic Imaging 2.4 (1993): 326-336.
"""
function mgradient(img; dims=coords_spatial(img), r=nothing, mode=:beucher)
    return mgradient(img, strel_box(img, dims; r); mode)
end
function mgradient(img::AbstractArray{T}, se; mode=:beucher) where {T}
    out = similar(img, maybe_floattype(T))
    buffer = _make_gradient_buffer(img, mode)
    return mgradient!(out, img, se, buffer; mode)
end

"""
    mgradient!(out, img, buffer; [dims], [r], [mode])
    mgradient!(out, img, se, buffer; [mode])

The in-place version of [`mgradient`](@ref) with input image `img` and output image `out`.

The `buffer` array is required for `:beucher` mode. For `:internal` and `:external` modes,
`buffer` is not needed and can be `nothing`.
"""
function mgradient!(out, img, buffer; dims=coords_spatial(img), r=nothing, mode=:beucher)
    return mgradient!(out, img, strel_box(img, dims; r), buffer; mode)
end

function mgradient!(out, img, se, @nospecialize(buffer); mode=:beucher)
    require_symmetric_strel(se)
    @debug "calculate mgradient using $mode mode"
    if mode == :beucher
        isnothing(buffer) && throw(ArgumentError("buffer array is required for mode :beucher"))
        _beucher_gradient!(out, img, se, buffer)
    elseif mode == :internal
        _internal_gradient!(out, img, se)
    elseif mode == :external
        _external_gradient!(out, img, se)
    end
    return out
end
mgradient!(out::AbstractArray{<:Color3}, img, se, buffer; kwargs...) = throw(ArgumentError("color image is not supported"))

_make_gradient_buffer(img, mode) = mode == :beucher ? similar(img, maybe_floattype(eltype(img))) : nothing

function _internal_gradient!(out, img, se)
    erode!(out, img, se)
    @. out = img - out
    return out
end

function _external_gradient!(out, img, se)
    dilate!(out, img, se)
    @. out = out - img
    return out
end

function _beucher_gradient!(out, img, se, buffer)
    dilate!(out, img, se)
    erode!(buffer, img, se)
    @. out = out - buffer
    return out
end
