# ---
# title: "`tophat`"
# id: op_tophat
# cover: assets/tophat.gif
# ---

# The (white) tophat operator is defined as `img - opening(img)`. Intuitively, this filter
# can be used to extract small white elements and details from an image.

using ImageMorphology
using TestImages
using ImageBase
using ImageShow

#-

img = restrict(testimage_dip3e("fig1005")) # wirebond mask

out1 = tophat(img) # default: all spatial dimensions, r=1, a box-shape SE
out2 = tophat(img; dims=(2,)) # only apply to the second dimension
out3 = tophat(img; r=5) # half-size r=5
out4 = tophat(img, strel_diamond((5, 5))) # generic SE -- this version doesn't accept `r` and `dims`

mosaic(img, out1, out2, out3, out4; nrow=1)


# ## See also

# To extract black small details, use [`bothat`](@ref op_bothat).

bothat(img) == tophat(complement.(img))

# For a comprehensive and more accurate documentation, please check the [`tophat`](@ref)
# reference page.

# save cover image #src
using ImageMagick #src
mkpath("assets") #src
outs = [tophat(img; r) for r in 1:5] #src
ImageMagick.save("assets/tophat.gif", cat(outs...; dims=3); fps=1) #src
