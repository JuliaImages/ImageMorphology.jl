"""
    filled_img = imfill(img::AbstractArray{Bool}, interval; dims=coords_spatial(img))
    filled_img = imfill(img::AbstractArray{Bool}, interval, connectivity)

Connected components of an image is found using flood-fill algorithm and returns a copy of
the original image after filling objects that falls in the range of interval.
For filling objects, represent the holes (part to be filled) with `true` in your array.

Parameters:
-  img            = Input image (Boolean array type)
-  interval       = objects of size (# of voxels) in this range will be filled with `false`
-  connectivity   = a Boolean-valued connectivity pattern, see [`label_components`](@ref).

# Examples

```jldoctest; setup=:(using ImageMorphology)
julia> img = Bool[0 0 1 1 0 0;
                  0 1 0 1 1 0;
                  0 0 1 1 0 0]
3×6 Matrix{Bool}:
 0  0  1  1  0  0
 0  1  0  1  1  0
 0  0  1  1  0  0

julia> imfill(.!(img), 0:3)
3×6 BitMatrix:
 1  1  0  0  1  1
 1  0  0  0  0  1
 1  1  0  0  1  1

julia> .!(ans)
3×6 BitMatrix:
 0  0  1  1  0  0
 0  1  1  1  1  0
 0  0  1  1  0  0
```
"""
function imfill(img::AbstractArray{Bool}, interval::Tuple{Real,Real}, connectivity::AbstractArray{Bool})
    return _imfill(img, interval, label_components(img, connectivity))
end
function imfill(img::AbstractArray{Bool}, interval::Tuple{Real,Real}; dims=coords_spatial(img))
    return _imfill(img, interval, label_components(img; dims=dims))
end
function imfill(img::AbstractArray{Bool}, interval, args...; kwargs...)
    return imfill(img, (minimum(interval)::Real, maximum(interval)::Real), args...; kwargs...)
end

function _imfill(img::AbstractArray{Bool}, interval::Tuple{Real,Real}, labels)
    if interval[1] > interval[2] || interval[1] < 0 || interval[2] < 0
        msg = "Interval must be non-negative and in format (min_range,max_range)"
        throw(DomainError(interval, msg))
    end
    count_labels = component_lengths(labels)

    new_img = similar(img)
    for ind in eachindex(img)
        if img[ind] == true && interval[1] <= count_labels[labels[ind]] <= interval[2]
            new_img[ind] = false
        else
            new_img[ind] = img[ind]
        end
    end

    return new_img
end
