# ---
# title: "`extreme_filter`"
# id: op_extreme_filter
# cover: assets/extreme_filter.png
# ---

# The `extreme_filter` function is the core operation in ImageMorphology. Many other
# morphological operations such as [`dilate`](@ref) and [`erode`](@ref) are direct usages of
# it. The cover image shows a pixel jitter using random select function.

# This page gives an overview of the `extreme_filter` function as a tutorial. It might not
# exactly show every details. For a comprehensive and more accurate documentation, please
# check the [`extreme_filter`](@ref) page.

using ImageMorphology
using ImageFiltering
using TestImages
using ImageBase
using ImageShow
using Random #hide
Random.seed!(1234); #hide

# ## Definition
#
# Below is the basic interface for this function
#
# ```julia
# extreme_filter(f, A, se)
# ```
#
# For each pixel `p` in A, a select function `f(x, y)` will be applied iteratively to
# its neighborhood. The neighborhood is defined by the structuring element `se`. In plain
# Julia codes, this is as simple as:
#
# ```julia
# # For illustration purposes only; you should use `extreme_filter` in practice
# function simple_extreme_filter(f, A, se)
#     se = strel(CartesianIndex, se)
#     R = CartesianIndices(A)
#     for p in R
#         x = A[p]
#         for o in se
#             q = p + o
#             q in R || continue
#             x = f(x, A[q])
#         end
#         A[p] = x
#     end
# end
# ```
#
# !!! tip "structuring element"
#     Structuring element defines the shape of each pixel's neighborhood. If you haven't
#     used [`strel`](@ref) before or don't know what `se` means. It is recommended to read
#     the concept introduction page "[Structuring element](@ref concept_se)" first.

# The name "extreme" comes from the fact that `max` and `min` are typical examples of the
# select function `f`.

# ## Examples

# To illustrate how it works, let's load the wirebond mask image and apply the extreme
# filter to it.

img = restrict(testimage_dip3e("fig1005")) # wirebond mask

#-

se = strel_diamond((3, 3)) # a dimond shape structuring element

# We use `max` as the select function -- applying `max` gives the [dilation](@ref op_dilate) operation:

extreme_filter(max, img, se)

# with the SE size controls the degree of dilation:

imgs = [extreme_filter(max, img, strel_diamond((n, n))) for n in 3:2:13]
mosaic(imgs; nrow=1)

# and `dims` controls the dimension to be dilated:

img_d1 = extreme_filter(max, img, strel_diamond((13, 13), (1,)))
img_d2 = extreme_filter(max, img, strel_diamond((13, 13), (2,)))
img_d12 = extreme_filter(max, img, strel_diamond((13, 13))) # (1, 2)
mosaic(img, img_d1, img_d2, img_d12; nrow=1)

# The `extreme_filter` can be used to build many operations, including the basic dilation,
# erosion, and some other operators. Here we just list a few for inspiration:

se = strel_diamond((7, 7))
img_1 = extreme_filter(max, img, se) # dilate
img_2 = extreme_filter(min, img, se) # erode
img_3 = extreme_filter((x, y) -> rand() > 0.5 ? x : y, img, se) # jitter
## or pipe the filters
img_4 = extreme_filter(max, extreme_filter(min, img), se) # opening
img_5 = extreme_filter(min, extreme_filter(max, img), se) # closing
mosaic(img, img_1, img_2, img_3, img_4, img_5; nrow=1)


# Many of the above operations are already implemented in ImageMorphology.

## the `mapwindow` function

# There is a very similar function to `extreme_filter` called `mapwindow(f, img, window)`
# from [ImageFiltering.jl](https://github.com/JuliaImages/ImageFiltering.jl). The main
# difference is that `f` in `mapwindow` applies to the entire neighborhood instead of two
# points in `extreme_filter`.

## currently, mapwindow only supports the box-shaped neighborhood
img1 = extreme_filter(max, img, strel_box((7, 7)))
img2 = mapwindow(maximum, img, (7, 7))
img1 == img2

# save cover image #src
using FileIO #src
mkpath("assets") #src
img = Gray.(restrict(testimage("blobs"))) #src
save("assets/extreme_filter.png", extreme_filter((x, y) -> rand() > 0.5 ? x : y, img; r=1)) #src
