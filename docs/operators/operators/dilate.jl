# ---
# title: "`dilate`"
# id: op_dilate
# cover: assets/dilate.gif
# ---

# The dilation operator `dilate` is essentially a max filter. This is a basic term in
# mathematical morphology -- many operations are built on top of `dilate` and the
# [`erode`](@ref) operator.

using ImageMorphology
using TestImages
using ImageBase
using ImageShow

img = restrict(testimage_dip3e("fig1005")) # wirebond mask

# For each pixel, dilation is the maximum of the pixels in the neighborhood. In mathematics,
# `dilate(A, Ω)` is defined as ``\delta_A[p] = \sup\{A[p+o] | o \in \Omega\}`` where `Ω` is
# the [structuring element](@ref concept_se).

# The simplest usage is `dilate(img; [dims], [r])`, where `dims` and `r` controls the
# neighborhood shape.

out1 = dilate(img) # default: all spatial dimensions, r=1, a box-shape SE
out2 = dilate(img; dims=(2,)) # only apply to the second dimension
out3 = dilate(img; r=5) # half-size r=5
mosaic(out1, out2, out3; nrow=1)

# It uses the [`strel_box`](@ref) to create a box-shaped structuring element. You can also
# provide a custom SE via the `dilate(img, se)` interface.

out1 = dilate(img, strel_box((3, 3))) # default se for`dilate(img)`
se = centered(Bool[1 1 1; 1 1 0; 0 0 0]) # match top-left region
out2 = dilate(img, se)
mosaic(out1, out2; nrow=1)

# An in-place version `dilate!` is also provided, for instance

out1 = similar(img)
dilate!(out1, img)

out2 = similar(img)
dilate!(out2, img, strel_diamond((3, 3)))
#md nothing #hide

# ## See also

# [`erode`](@ref op_erode) is the dual operator of `dilate` in the following sense:

complement.(dilate(img)) == erode(complement.(img))

# For bool arrays and [symmetric SEs](@ref concept_symmetric), dilation becomes equivalent
# to the minkowski addition on sets: ``A \oplus B = \{ a+b | a \in A, b \in B \}``.

# For a comprehensive and more accurate documentation, please check the [`dilate`](@ref)
# reference page.

# save cover image #src
using ImageMagick #src
mkpath("assets") #src
outs = [dilate(img; r) for r in 1:2:9] #src
ImageMagick.save("assets/dilate.gif", cat(outs...; dims=3); fps=1) #src
