using ImageMorphology
using TestImages
using ImageBase
using ImageShow

img = restrict(testimage_dip3e("fig1005")) # wirebond mask

out1 = tophat(img) # default: all spatial dimensions, r=1, a box-shape SE
out2 = tophat(img; dims=(2,)) # only apply to the second dimension
out3 = tophat(img; r=5) # half-size r=5
out4 = tophat(img, strel_diamond((5, 5))) # generic SE -- this version doesn't accept `r` and `dims`

mosaic(img, out1, out2, out3, out4; nrow=1)

bothat(img) == tophat(complement.(img))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

