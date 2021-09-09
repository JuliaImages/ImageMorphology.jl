"""
 Image filling Algorithm

 ```
 filled_img = imfill(img, interval)
 filled_img = imfill(img, interval, value)
 filled_img = imfill(img, interval, value, connectivity)
 ```

 Connected components of an image is found using flood-fill algorithm and returns a copy of
 the original image after filling objects that falls in the range of interval.
 For filling objects, represent the holes(part to be filled) with `true` in your array.


 Parameters:

  -  img            = Input image (Boolean array type)
  -  interval       = objects of size in this range will be filled with `false`
  -  connectivity   = connectivity takes the same values as in label_components (Default value is 1:ndims(img))

 """

function imfill(img::AbstractArray{Bool}, interval::Tuple{Real,Real}, connectivity::Union{Dims, AbstractVector{Int}, BitArray}=1:ndims(img))

    if interval[1] > interval[2] || interval[1] < 0 || interval[2] < 0
        throw(DomainError(interval,"Interval must be non-negative and in format (min_range,max_range)"))
    end

    labels = label_components(img,connectivity)

    count_labels = component_lengths(labels)

    new_img = similar(img)
    for ind in eachindex(img)
        if img[ind] == true && interval[1] <= count_labels[labels[ind]+1] <= interval[2]
            new_img[ind] = false
        else
            new_img[ind] = img[ind]
        end
    end

    return new_img
 end
