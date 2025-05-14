# ---
# title: "`mreconstruct`"
# id: op_mreconstruct
# cover: assets/mreconstruct.png
# ---

# The morphological reconstruction operator is to repeatedly apply particular operator until
# the stability, i.e,. output unchanged. The most widely used ones are reconstruction by
# dilation and reconstruction by erosion.

using ImageMorphology
using TestImages
using ImageBase
using ImageShow
using IndirectArrays, FileIO #hide

img = Gray.(testimage("blob"))
nothing #hide

# Unlike other morphological operators such as `erode`, `mreconstruct` requires two image
# inputs. One is the marker image, and the other is the mask image. Usually, `mask` is the
# image we study on and the `marker` image is what we want to reconstruct from. A suitable
# marker image can be determined using[^1]:

# - knowledge about the expected result
# - known facts about the image or the physics of the object it represents
# - some transformations of the mask image itself
# - other image data, if available (i.e., multispectral and multitemporal images)
# - interaction with the user (i.e., markers are manually defined)

# Take reconstruction by dilation as an example, we build the marker image from the mask
# image, then repeatedly apply `marker = min.(dilate(marker), mask)`.

mask = img
marker = mask .- 0.3
mosaic(marker, mask; nrow=1)

# The following function is a naive implementation of this idea:

# ```julia
# # For illustration purposes only; you should use `mreconstruct` in practice
# function my_reconstruct(marker, mask)
#     # reconstruction by dilation
#     # For reconstruction by erosion, one can change `min` to `max` and `dilate` to `erode`
#     out = min.(mask, marker)
#     has_changed = true
#     while has_changed
#         old = out
#         out = min.(dilate(out), mask)
#         has_changed = !(old == out)
#     end
#     return out
# end
# ```

# The following GIF shows the intermediate result `out` of each iteration. The left side
# shows its actual gray value, and the right side shows the colored value to better
# visualize the progress.

function my_reconstruct_outs(marker, mask) #hide
    out = min.(mask, marker) #hide
    outs = [out] #hide
    has_changed = true #hide
    n = 1 #hide
    while has_changed #hide
        old = out #hide
        out = min.(dilate(out), mask) #hide
        n % 4 == 1 && push!(outs, out) #hide
        has_changed = !(old == out) #hide
        n = n + 1 #hide
    end #hide
    return outs #hide
end #hide
function with_color(img, orders, colormap) #hide
    return IndirectArray([orders[img[i]] for i in CartesianIndices(img)], colormap) #hide
end #hide
outs = my_reconstruct_outs(marker, mask) #hide
orders = Dict(x => i for (i, x) in enumerate(sort(mapreduce(unique, union, outs)))) #hide
colormap = ImageCore.Colors.distinguishable_colors(length(orders)) #hide
f = scaleminmax(RGB, 0, 1) #hide
save("assets/mreconstruct.gif",  #hide
     f.(cat([mosaic(frame, with_color(frame, orders, colormap); nrow=1) #hide
             for frame in (@view outs[1:1:end])]..., dims=3)); #hide
     fps=5) #hide
# ![alternative text](assets/mreconstruct.gif)

# This is how reconstruction happends and why it is called so -- it reconstructs the
# `marker` image with a reference image `mask`, this helps us investigate certain properties
# of the `mask` image. For instance, if we want to investigate the blobs only (dark region),
# we would want to apply reconstruction by dilation to it.

# The `mreconstruct(op, marker, mask)` function currently supports only `dilate` and `erode`
# as `op`, indicating reconstruction by dilation and reconstruction by erosion,
# respectively. There are also aliases [`underbuild`](@ref op_underbuild) and
# [`overbuild`](@ref op_overbuild).

## reconstruction by dilation -- `underbuild`
out1 = mreconstruct(dilate, img .- 0.3, img)
## reconstruction by erosion -- `overbuild`
out2 = mreconstruct(erode, img .+ 0.3, img)
mosaic(mask, out1, out2; nrow=1)

# The default structuring element it uses is the box-shaped SE. Like the other operators,
# you can pass `dims` keyword or positional argument `se` to `mreconstruct`. For instance,

mreconstruct(dilate, img .- 0.3, img; dims=1) # reconstruction along colomn
mreconstruct(dilate, img .- 0.3, img, strel_diamond(img)) # diamond shape SE
nothing #hide

# However, for generic structuring element, the half-size (or, radius) should be at most `1`
# for each dimension. This is because the algorithm is a growing process until stability,
# thus SE with half-size larger than `1` doesn't make a difference in practice.

# save cover image #src
using FileIO #src
img = Gray.(restrict(testimage("blobs"))) #src
save("assets/mreconstruct.png", mreconstruct(dilate, img .- 0.3, img)) #src

# [^1]: P. Soille, Morphological Image Analysis (Chapter 6.2.3). Berlin, Heidelberg: Springer Berlin Heidelberg, 2004. doi: 10.1007/978-3-662-05088-0.
