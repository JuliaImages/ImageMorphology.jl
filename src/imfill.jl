"""
 Image filling Algorithm

 ```
 filled_img = imfill(img, min_size)
 filled_img = imfill(img, min_size, threshold)
 filled_img = imfill(img, min_size, threshold, value)
 filled_img = imfill(img, min_size, threshold, value, connectivity)
 ```

 Connected components of an image is found using flood-fill algorithm and returns a copy of
 the original image after filling objects that are smaller than the specified min_size.


 Parameters:

  -  img            = Input image (Float valued)
  -  min_size       = Minimum size of objects that are not filled.
  -  threshold      = Pixels below this value will be ignored while finding connected components (Default value is 0)
  -  value          = Value to be given to the objects that are filled (Default value is 1)
  -  connectivity   = connectivity takes the same values as in label_components (Default value is 1:ndims(img))

 """

  function imfill(img::AbstractArray{float}, min_size::Int, threshold::AbstractFloat=0.0, value::Int=1, connectivity::Union{Dims, AbstractVector{Int}, BitArray}=1:ndims(img))
     if min_size == 0 || min_size == 1
         return img
     end

     to_bool(x::AbstractFloat) = x<= threshold ? false : x>0 ? true : throw(InexactError())
     bool_img = to_bool.(img)
     inverted = .!(bool_img)

     labels = label_components(inverted,connectivity)

     number_labels = unique(labels)
     count_labels = Dict([(i,count(x->x==i,labels)) for i in number_labels])

     new_img = similar(img)
     for itr in eachindex(img)
         if img[itr] == 0
             if count_labels[labels[itr]] < min_size
                 new_img[itr] = value
             else
                 new_img[itr] = 0
             end
         else
             new_img[itr] = img[itr]
         end
     end

      return new_img
 end
