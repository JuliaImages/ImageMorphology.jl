# ---
# title: "`underbuild`"
# id: op_underbuild
# cover: assets/underbuild.png
# ---

# `underbuild` is reconstruction by dilation, it is an alias for `mreconstruct` when `op=dilate`.

using ImageMorphology
using TestImages
using ImageBase
using ImageShow

img = restrict(Gray.(testimage("blob")))
nothing # hide

# `underbuild` accepts two images as inputs: `marker` and `mask`.

out = underbuild(img .- 0.2, img)
mosaic(img, out; nrow=1)

# ## See also

# [`underbuild`](@ref op_underbuild) is the dual operator of `overbuild` in the following sense:

marker, mask = rand(32, 32), rand(32, 32)
complement.(underbuild(marker, mask)) == underbuild(complement.(marker), complement.(mask))

# For more details, please refer to [`mreconstruct`](@ref op_mreconstruct).

# save cover image #src
using FileIO #src
mkpath("assets") #src
img = Gray.(restrict(testimage("blobs"))) #src
save("assets/underbuild.png", underbuild(img .- 0.2, img)) #src
