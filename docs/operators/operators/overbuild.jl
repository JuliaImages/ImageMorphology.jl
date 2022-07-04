# ---
# title: "`overbuild`"
# id: op_overbuild
# cover: assets/overbuild.png
# ---

# `overbuild` is reconstruction by erosion, it is an alias for `mreconstruct` when `op=erode`.

using ImageMorphology
using TestImages
using ImageBase
using ImageShow

img = restrict(Gray.(testimage("blob")))
nothing # hide

# `overbuild` accepts two images as inputs: `marker` and `mask`.

out = overbuild(img .- 0.2, img)
mosaic(img, out; nrow=1)

# ## See also

# [`overbuild`](@ref op_overbuild) is the dual operator of `underbuild` in the following sense:

marker, mask = rand(32, 32), rand(32, 32)
complement.(overbuild(marker, mask)) == overbuild(complement.(marker), complement.(mask))

# For more details, please refer to [`mreconstruct`](@ref op_mreconstruct).

# save cover image #src
using FileIO #src
mkpath("assets") #src
img = Gray.(restrict(testimage("blobs"))) #src
save("assets/overbuild.png", overbuild(img .- 0.2, img)) #src
