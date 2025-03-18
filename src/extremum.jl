using ImageCore.MappedArrays

"""
    hmaxima(img, h; [dims])
    hmaxima(img, h; se)

Compute regional maxima of height h.

For grayscale image `img`, the hmaxima transformation suppresses all regional maxima whose
height is below a given threshold level `h`.

The `dims` keyword is used to specify the dimension to process by constructing the box shape
structuring element [`strel_box(img; dims)`](@ref strel_box). For generic structuring
element, the half-size is expected to be either `0` or `1` along each dimension.

# Examples

```jldoctest; setup=:(using ImageMorphology)
julia> img = UInt8[
    3 3 3 3 3 3 3 3 3 3
    3 4 4 4 3 3 4 4 4 3
    3 4 5 4 3 3 4 9 4 3
    3 4 4 4 3 3 4 4 4 3
    3 3 3 3 3 3 3 3 3 3
    ]
julia> out = hmaxima(reinterpret(N0f8, img), reinterpret(N0f8, UInt8(3)))
julia> ref_img = UInt8[
    3 3 3 3 3 3 3 3 3 3
    3 3 3 3 3 3 4 4 4 3
    3 3 3 3 3 3 4 6 4 3
    3 3 3 3 3 3 4 4 4 3
    3 3 3 3 3 3 3 3 3 3
    ]
julia> @test eltype(out) == N0f8
julia> @test out == reinterpret(N0f8, ref_img)
```
    
# See also

The inplace version of this function is `hmaxima!`.

# References

- [1] L. Vincent, “Morphological grayscale reconstruction in image analysis: applications
  and efficient algorithms,” IEEE Trans. on Image Process., vol. 2, no. 2, pp. 176–201, Apr.
  1993, doi: 10.1109/83.217222.
- [2] P. Soille, Morphological Image Analysis. Berlin, Heidelberg: Springer Berlin
  Heidelberg, 2004. doi: 10.1007/978-3-662-05088-0.
"""
function hmaxima(img, h, se)
    return _hmaxima!(similar(img), img, h, se)
end

function hmaxima(img, h; dims=coords_spatial(img))
    return _hmaxima!(similar(img), img, h, strel_box(img, dims))
end

function hmaxima!(out, img, h; dims=coords_spatial(img))
    return _hmaxima!(out, img, h, strel_box(img, dims))
end

function _hmaxima!(out, img, h, se)
    # tmp = map(x -> saturating_sub(x, h), img) 
    tmp = mappedarray(x -> saturating_sub(x, h), img)
    return mreconstruct!(dilate, out, reinterpret(eltype(img), tmp), img, se)
end

"""
    hminima(img, h; [dims])
    hminima(img, h; se)

Compute regional minima of depth h.

For grayscale image `img`, the hminima transformation suppresses all regional minima whose
depth is below a given threshold level `h`.

The `dims` keyword is used to specify the dimension to process by constructing the box shape
structuring element [`strel_box(img; dims)`](@ref strel_box). For generic structuring
element, the half-size is expected to be either `0` or `1` along each dimension.

# Examples

```jldoctest; setup=:(using ImageMorphology)
julia> img = UInt8[
        8 8 8 8 8 8 8 8 8 8
        8 7 7 7 8 8 7 7 7 8
        8 7 6 7 8 8 7 3 7 8
        8 7 7 7 8 8 7 7 7 8
        8 8 8 8 8 8 8 8 8 8
    ]
julia> out = hminima(reinterpret(N0f8, img), reinterpret(N0f8, UInt8(3)))
julia> ref_img = UInt8[
        8 8 8 8 8 8 8 8 8 8
        8 8 8 8 8 8 7 7 7 8
        8 8 8 8 8 8 7 6 7 8
        8 8 8 8 8 8 7 7 7 8
        8 8 8 8 8 8 8 8 8 8
    ]
julia> @test eltype(out) == N0f8
julia> @test out == reinterpret(N0f8, ref_img)

# See also

The inplace version of this function is `hminima!`.

# References

- [1] L. Vincent, “Morphological grayscale reconstruction in image analysis: applications
  and efficient algorithms,” IEEE Trans. on Image Process., vol. 2, no. 2, pp. 176–201, Apr.
  1993, doi: 10.1109/83.217222.
- [2] P. Soille, Morphological Image Analysis. Berlin, Heidelberg: Springer Berlin
  Heidelberg, 2004. doi: 10.1007/978-3-662-05088-0.
"""
function hminima(img, h, se)
    return _hminima!(similar(img), img, h, se)
end

function hminima(img, h; dims=coords_spatial(img))
    return _hminima!(similar(img), img, h, strel_box(img, dims))
end

function hminima!(out, img, h; dims=coords_spatial(img))
    return _hminima!(out, img, h, strel_box(img, dims))
end

function _hminima!(out, img, h, se)
    # tmp = map(x -> saturating_add(x, h), img)
    tmp = mappedarray(x -> saturating_add(x, h), img)
    return mreconstruct!(erode, out, reinterpret(eltype(img), tmp), img, se)
end

"""
    regional_maxima(img; [dims])
    regional_maxima(img; se)

Determines all regional maxima of the `image`.

For grayscale image `img`, A regional maximum is defined as the connected set of pixels that have the
same value, which is greater than the values of all pixels in direct neighborhood of the set.

The `dims` keyword is used to specify the dimension to process by constructing the box shape
structuring element [`strel_box(img; dims)`](@ref strel_box). For generic structuring
element, the half-size is expected to be either `0` or `1` along each dimension.

#Note
This implementation is faster than maxtree local_maxima approach if maxtree not precomputed

# See also

The inplace version of this function is `regional_maxima!`.

# References

- [1] L. Vincent, “Morphological grayscale reconstruction in image analysis: applications
    and efficient algorithms,” IEEE Trans. on Image Process., vol. 2, no. 2, pp. 176–201, Apr.
    1993, doi: 10.1109/83.217222.
- [2] P. Soille, Morphological Image Analysis. Berlin, Heidelberg: Springer Berlin
    Heidelberg, 2004. doi: 10.1007/978-3-662-05088-0.
"""
function regional_maxima(img, se)
    return _regional_maxima!(similar(img, Bool), img, se)
end

function regional_maxima(img; dims=coords_spatial(img))
    return _regional_maxima!(similar(img, Bool), img, strel_box(img, dims))
end

function regional_maxima!(out, img; dims=coords_spatial(img))
    return _regional_maxima!(out, img, strel_box(img, dims))
end

function _regional_maxima!(out, img, se)
    #print(reinterpret(eltype(img), (Unsigned)1))
    h = reinterpret(eltype(img), (UInt8)(1))
    out = ((img - hmaxima(img, h, se))) .> 0.0
    return out
end


"""
    regional_minima(img; [dims])
    regional_minima(img; se)

Deetermines all regional minimum of the `image`.

For grayscale image `img`, A regional minimum is defined as the connected set of pixels that have the
same value, which is lower than the values of all pixels in direct neighborhood of the set.

The `dims` keyword is used to specify the dimension to process by constructing the box shape
structuring element [`strel_box(img; dims)`](@ref strel_box). For generic structuring
element, the half-size is expected to be either `0` or `1` along each dimension.

#Note
This implementation is faster than maxtree local_minima approach if maxtree not precomputed

# See also

The inplace version of this function is `regional_minima!`.

# References

- [1] L. Vincent, “Morphological grayscale reconstruction in image analysis: applications
    and efficient algorithms,” IEEE Trans. on Image Process., vol. 2, no. 2, pp. 176–201, Apr.
    1993, doi: 10.1109/83.217222.
- [2] P. Soille, Morphological Image Analysis. Berlin, Heidelberg: Springer Berlin
    Heidelberg, 2004. doi: 10.1007/978-3-662-05088-0.
"""
function regional_minima(img, se)
    return _regional_minima!(similar(img, Bool), img, se)
end

function regional_minima(img; dims=coords_spatial(img))
    return _regional_minima!(similar(img, Bool), img, strel_box(img, dims))
end

function regional_minima!(out, img; dims=coords_spatial(img))
    return _regional_minima!(out, img, strel_box(img, dims))
end

function _regional_minima!(out, img, se)
    h = reinterpret(eltype(img), (UInt8)(1))
    out = ((hminima(img, h, se) - img)) .> 0.0
    return out
end
