using ImageMorphology
using TestImages
using ImageBase
using ImageShow

img = restrict(testimage_dip3e("fig0940")) # rice
img01 = @. Gray(img > 0.5)
mosaic(img, img01; nrow=1)

out1 = closing(img) # default: all spatial dimensions, r=1, a box-shape SE
out2 = closing(img; dims=(2,)) # only apply to the second dimension
out3 = closing(img; r=5) # half-size r=5

# also to the binary version
out1_01 = closing(img01)
out2_01 = closing(img01; dims=(2,))
out3_01 = closing(img01; r=5)

mosaic(out1, out2, out3, out1_01, out2_01, out3_01; nrow=2, rowmajor=true)

complement.(closing(img)) == closing(complement.(img))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

