"""
 Image filling Algorithm

 ```
 filled_img = imfill(img, min_size)
 filled_img = imfill(img, min_size, value)
 filled_img = imfill(img, min_size, value, connectivity)
 ```

 Connected components of an image is found using flood-fill algorithm and returns a copy of
 the original image after filling objects that are smaller than the specified min_size.


 Parameters:

  -  img            = Input image (Boolean array type)
  -  min_size       = Minimum size of objects that are not filled.
  -  value          = Value to be given to the objects that are filled (Default value is 1)
  -  connectivity   = connectivity takes the same values as in label_components (Default value is 1:ndims(img))

 """

  function imfill(img::BitArray, min_size::Int, threshold::AbstractFloat=0.0, value::Int=1, connectivity::Union{Dims, AbstractVector{Int}, BitArray}=1:ndims(img))
     if min_size == 0 || min_size == 1
         return img
     end

     labels = label_components(img,connectivity)

     number_labels = unique(labels)
     count_labels = Dict([(i,count(x->x==i,labels)) for i in number_labels])

     new_img = similar(img)
     for ind in eachindex(img)
         if img[ind] == 0
             if count_labels[labels[ind]] < min_size
                 new_img[ind] = value
             else
                 new_img[ind] = 0
             end
         else
             new_img[ind] = img[ind]
         end
     end

      return new_img
 end
