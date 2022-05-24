# ---
# title: "`mgradient`"
# id: op_mgradient
# cover: assets/mgradient.png
# ---

# There are three commonly used morphological gradient definitions: the beucher gradident
# `dilate(img) - erode(img)`, the external half-gradient `dilate(img) - img`, and the
# internal half-gradident `img - erode(img)`.

using ImageMorphology
using TestImages
using ImageBase
using ImageShow

img = restrict(testimage_dip3e("fig0940")) # rice
img01 = @. Gray(img > 0.5)
mosaic(img, img01; nrow=1)

# ## Beucher gradient

# The default mode is beucher gradient defined as `dilate(img) - erode(img)`.

mosaic(mgradient(img), mgradient(img01); nrow=1)

# There are two nice properties of the beucher gradient:
#
# - it is always positive if the structuring element is symmetric
# - it is a self-complementary operator -- `mgradient(img) == mgradient(complement.(img))`

all(mgradient(rand(32, 32)) .> 0) # always positive

#-

A = rand(32, 32)
mgradient(A) == mgradient(1 .- A) # self-complementary

# ## Half gradient

# If we look closely, we can see that the beucher gradient boundary is always two-pixel wide
# -- one pixel per each side.

A = falses(7, 7)
A[3:5, 3:5] .= true
A

#-

Int.(mgradient(A))

# This is because `dilate` extends the boundary for one pixel, while `erode` shrinks for one
# pixel. This observation introduces so-called _half-gradients_.
#
# The _external half-gradient_, also known as _half-gradient by dilation_, is defined as
# `dilate(img) - img`:

Int.(mgradient(A; mode=:external)) # external boundary

# The _internal half-gradient_, also known as _half-gradient by erosion_, is defined as
# `img - erode(img)`:

Int.(mgradient(A; mode=:internal)) # internal boundary

# Note that the result of internal gradient is inside the original image, and that of
# external gradient is outside the original image. Also, external gradient and internal
# are complementary to each other:

A = rand(32, 32)
mgradient(A; mode=:external) == mgradient(1 .- A; mode=:internal)

# ## Thick gradient

# The gradient thickness is controlled by the SE size, or the `r` parameter. When `r>1`, it
# is usually called _thick gradient_.

mosaic([mgradient(img01; r) for r in 1:3]; nrow=1)

# ## See also

# There are many useful gradient definitions in image processing. For instance, ImageBase.jl
# provides the finite-difference version of gradients `fdiff` and `fgradient`.
# ImageFiltering.jl provides `sobel`, `prewitt` and other gradient filters.

fout = abs.(ImageBase.FiniteDiff.fdiff(img; dims=1))
mout = mgradient(img; dims=1)
mosaic(fout, mout; nrow=1)


# For a comprehensive and more accurate documentation, please check the [`mgradient`](@ref)
# reference page.

# save cover image #src
using FileIO #src
mkpath("assets") #src
img = Gray.(restrict(testimage("blobs"))) #src
save("assets/mgradient.png", mgradient(img; r=2)) #src
