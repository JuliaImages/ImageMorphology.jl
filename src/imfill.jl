"""
 Image filling Algorithm

 ```
 filled_img = imfill(img, interval)
 filled_img = imfill(img, interval, value)
 filled_img = imfill(img, interval, value, connectivity)
 ```

 Connected components of an image is found using flood-fill algorithm and returns a copy of
 the original image after filling objects that falls in the range of interval.


 Parameters:

  -  img            = Input image (Boolean array type)
  -  interval       = objects of size in this range will be filled with False (default value is (0,64) )
  -  connectivity   = connectivity takes the same values as in label_components (Default value is 1:ndims(img))

 """

  function imfill(img::AbstractArray{Bool}, interval::Tuple{Real,Real}=(0,64), connectivity::Union{Dims, AbstractVector{Int}, BitArray}=1:ndims(img))

     labels = label_components(img,connectivity)

     count_labels = Dict([(i,count(x->x==i,labels)) for i in unique(labels)])

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
