"""
```
cleared_img = clearborder(img)
cleared_img = clearborder(img, width)
cleared_img = clearborder(img, width, background)
```
Returns a copy of the original image after clearing objects connected to the border of the image.
Parameters:
 -  img          = Binary/Grayscale input image
 -  width        = Width of the border examined (Default value is 1)
 -  background   = Value to be given to pixels/elements that are cleared (Default value is 0)
"""
function clearborder(img::AbstractArray, width::Integer=1, background::Integer=0)

    for i in size(img)
        if (width > i)
            throw(ArgumentError("Border width must not be greater than size of the image."))
        end
    end

    labels = label_components(img, strel_box(img))
    number = maximum(labels) + 1

    dimensions = size(img)
    outerrange = CartesianIndices(map(i -> 1:i, dimensions))
    innerrange = CartesianIndices(map(i -> (1 + width):(i - width), dimensions))

    border_labels = Set{Int}()
    for i in EdgeIterator(outerrange, innerrange)
        push!(border_labels, labels[i])
    end

    new_img = similar(img)
    for itr in eachindex(labels)
        if labels[itr] in border_labels
            new_img[itr] = background
        else
            new_img[itr] = img[itr]
        end
    end

    return new_img

end
