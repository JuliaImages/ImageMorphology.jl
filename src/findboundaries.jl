"""
```
bordered_img = findboundries(img)
bordered_img = findboundries(img, background)
```
Returns a copy of the original image where boundries between labeled regions are true.
Parameters:
 -  img             = Binary/Grayscale input image
 -  background      = Value to be given to pixels/elements that are cleared (Default value is 0)
"""

function findboundaries(img::AbstractArray, background::Integer=0)

    connectivity = ntuple(i -> 3, ndims(img))
    labels = label_components(img,trues(connectivity))

    boundaries = dilate(img) .!= erode(img)
    foreground_image = (labels .!= background)
    boundaries = boundaries .& foreground_image
    
    return boundaries
end
