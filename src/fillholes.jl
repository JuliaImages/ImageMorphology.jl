"""
    fillhole(img; [dims])
    fillhole(img; se)

Fill the holes in image 'img'. Could be binary or grascale

The `dims` keyword is used to specify the dimension to process by constructing the box shape
structuring element [`strel_box(img; dims)`](@ref strel_box). For generic structuring
element, the half-size is expected to be either `0` or `1` along each dimension.

The output has the same type as input image
"""

function fillhole(img; dims=coords_spatial(img))
    return fillhole(img, strel_box(img, dims))
end

function fillhole(img, se)
    return fillhole!(similar(img), img, se)
end

function fillhole!(out, img; dims=coords_spatial(img))
    return fillhole!(out, img, strel_box(img, dims))
end

function fillhole!(out, img, se)
    return _fillhole!(out, img, se)
end

function _fillhole!(out, img, se)
    N = ndims(img)

    axes(out) == axes(img) || throw(DimensionMismatch("images should have the same axes"))

    se_size = strel_size(se)
    if length(se_size) != N
        msg = "the input structuring element is not for $N dimensional array, instead it is for $(length(se_size)) dimensional array"
        throw(DimensionMismatch(msg))
    end
    if !all(x -> in(x, (1, 3)), strel_size(se))
        msg = "structuring element with half-size larger than 1 is invalid"
        throw(DimensionMismatch(msg))
    end

    tmp = similar(img)

    # fill marker image with max
    fill!(tmp, typemax(eltype(img)))
    # fill borders with 0
    dimensions = size(tmp)
    outerrange = CartesianIndices(map(i -> 1:i, dimensions))
    innerrange = CartesianIndices(map(i -> (1 + 1):(i - 1), dimensions))
    for i in EdgeIterator(outerrange, innerrange)
        tmp[i] = 0
    end

    return mreconstruct!(erode, out, tmp, img, se)
end