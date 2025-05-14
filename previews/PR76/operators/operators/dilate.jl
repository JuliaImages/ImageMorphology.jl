using ImageMorphology
using TestImages
using ImageBase
using ImageShow

img = restrict(testimage_dip3e("fig1005")) # wirebond mask

out1 = dilate(img) # default: all spatial dimensions, r=1, a box-shape SE
out2 = dilate(img; dims=(2,)) # only apply to the second dimension
out3 = dilate(img; r=5) # half-size r=5
mosaic(out1, out2, out3; nrow=1)

out1 = dilate(img, strel_box((3, 3))) # default se for`dilate(img)`
se = centered(Bool[1 1 1; 1 1 0; 0 0 0]) # match top-left region
out2 = dilate(img, se)
mosaic(out1, out2; nrow=1)

out1 = similar(img)
dilate!(out1, img)

out2 = similar(img)
dilate!(out2, img, strel_diamond((3, 3)))

complement.(dilate(img)) == erode(complement.(img))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

