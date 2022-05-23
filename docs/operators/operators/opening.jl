# ---
# title: "`opening`"
# id: op_opening
# cover: assets/opening.gif
# ---

# `opening` operator is defined as `dilate(erode(img))`. Intuitively, opening operation
# fills the white holes in the image.

using ImageMorphology
using TestImages
using ImageBase
using ImageShow

img = restrict(testimage_dip3e("fig0940")) # rice
img01 = @. Gray(img > 0.5)
mosaic(img, img01; nrow=1)

#-

out1 = opening(img)
out2 = opening(img; dims=(2,))
out3 = opening(img; r=5)

out1_01 = opening(img01)
out2_01 = opening(img01; dims=(2,))
out3_01 = opening(img01; r=5)

mosaic(out1, out2, out3, out1_01, out2_01, out3_01; nrow=2, rowmajor=true)


# ## See also

# [`closing`](@ref op_closing) is the dual operator of `opening` in the following sense:

complement.(opening(img)) == closing(complement.(img))

# For a comprehensive and more accurate documentation, please check the [`opening`](@ref)
# reference page.

# save cover image #src
using ImageMagick #src
mkpath("assets") #src
outs = [opening(img01; r) for r in 1:5] #src
ImageMagick.save("assets/opening.gif", cat(outs...; dims=3); fps=1) #src
