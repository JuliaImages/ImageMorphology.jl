# ---
# title: "`erode`"
# id: op_erode
# cover: assets/erode.gif
# ---

# The dilation operator `erode` is essentially a min filter. This is a basic term in
# mathematical morphology -- many operations are built on top of `erode` and the
# [`dilate`](@ref) operator.

using ImageMorphology
using TestImages
using ImageBase
using ImageShow

img = restrict(testimage_dip3e("fig1005")) # wirebond mask

# For each pixel, dilation is the minimum of the pixels in the neighborhood. In mathematics,
# `erode(A, Ω)` is defined as ``\delta_A[p] = \inf\{A[p+o] | o \in \Omega\}`` where `Ω` is
# the [structuring element](@ref concept_se).

# The simplest usage is `erode(img; [dims], [r])`, where `dims` and `r` controls the
# neighborhood shape.

out1 = erode(img) # default: all spatial dimensions, r=1, a box-shape SE
out2 = erode(img; dims=(2,)) # only apply to the second dimension
out3 = erode(img; r=5) # half-size r=5
mosaic(out1, out2, out3; nrow=1)

# It uses the [`strel_box`](@ref) to create a box-shaped structuring element. You can also
# provide a custom SE via the `erode(img, se)` interface.

out1 = erode(img, strel_box((3, 3))) # default se for`erode(img)`
se = centered(Bool[1 1 1; 1 1 0; 0 0 0]) # match top-left region
out2 = erode(img, se)
mosaic(out1, out2; nrow=1)

# An in-place version `erode!` is also provided, for instance

out1 = similar(img)
erode!(out1, img)

out2 = similar(img)
erode!(out2, img, strel_diamond((3, 3)))
#md nothing #hide

# ## See also

# [`erode`](@ref op_erode) is the dual operator of `dilate` in the following sense:

complement.(erode(img)) == erode(complement.(img))

# For bool arrays and [symmetric SEs](@ref concept_symmetric), erosion becomes equivalent
# to the minkowski difference on sets: ``A \ominus B = \{ a-b | a \in A, b \in B \}``.

# For a comprehensive and more accurate documentation, please check the [`erode`](@ref)
# reference page.

# save cover image #src
using ImageMagick #src
mkpath("assets") #src
outs = [erode(img; r) for r in 1:2:9] #src
ImageMagick.save("assets/erode.gif", cat(outs...; dims=3); fps=1) #src
