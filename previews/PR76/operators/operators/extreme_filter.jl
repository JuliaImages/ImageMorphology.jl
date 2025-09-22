using ImageMorphology
using ImageFiltering
using TestImages
using ImageBase
using ImageShow
using Random #hide
Random.seed!(1234); #hide

img = restrict(testimage_dip3e("fig1005")) # wirebond mask

se = strel_diamond((3, 3)) # a dimond shape structuring element

extreme_filter(max, img, se)

imgs = [extreme_filter(max, img, strel_diamond((n, n))) for n in 3:2:13]
mosaic(imgs; nrow=1)

img_d1 = extreme_filter(max, img, strel_diamond((13, 13), (1,)))
img_d2 = extreme_filter(max, img, strel_diamond((13, 13), (2,)))
img_d12 = extreme_filter(max, img, strel_diamond((13, 13))) # (1, 2)
mosaic(img, img_d1, img_d2, img_d12; nrow=1)

se = strel_diamond((7, 7))
img_1 = extreme_filter(max, img, se) # dilate
img_2 = extreme_filter(min, img, se) # erode
img_3 = extreme_filter((x, y) -> rand() > 0.5 ? x : y, img, se) # jitter
# or pipe the filters
img_4 = extreme_filter(max, extreme_filter(min, img), se) # opening
img_5 = extreme_filter(min, extreme_filter(max, img), se) # closing
mosaic(img, img_1, img_2, img_3, img_4, img_5; nrow=1)

# the `mapwindow` function

# currently, mapwindow only supports the box-shaped neighborhood
img1 = extreme_filter(max, img, strel_box((7, 7)))
img2 = mapwindow(maximum, img, (7, 7))
img1 == img2

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

