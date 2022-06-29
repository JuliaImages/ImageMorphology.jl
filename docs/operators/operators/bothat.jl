# ---
# title: "`bothat`"
# id: op_bothat
# cover: assets/bothat.png
# ---

# The (black) tophat operator, also known as bottom hat, is defined as `closing(img) - img`.
# Intuitively, this filter can be used to extract small black elements and details from an
# image.

using ImageMorphology
using TestImages
using ImageBase
using ImageShow

#-

img = restrict(testimage_dip3e("fig1038")) # fingerprint

out1 = bothat(img) # default: all spatial dimensions, r=1, a box-shape SE
out2 = bothat(img; dims=(2,)) # only apply to the second dimension
out3 = bothat(img; r=5) # half-size r=5
out4 = bothat(img, strel_diamond((5, 5))) # generic SE -- this version doesn't accept `r` and `dims`

mosaic(img, out1, out2, out3, out4; nrow=1)

# ## See also

# To extract white small details, use the white [`tophat`](@ref op_tophat):

tophat(img) == bothat(complement.(img))

# For a comprehensive and more accurate documentation, please check the [`bothat`](@ref)
# reference page.

# save cover image #src
using FileIO #src
mkpath("assets") #src
img = Gray.(restrict(testimage("blobs"))) #src
save("assets/bothat.png", bothat(img; r=2)) #src
