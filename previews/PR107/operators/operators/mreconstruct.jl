using ImageMorphology
using TestImages
using ImageBase
using ImageShow
using IndirectArrays #hide

img = Gray.(testimage("blob"))
nothing #hide

mask = img
marker = mask .- 0.3
mosaic(marker, mask; nrow=1)

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
gif = ImageShow.gif([mosaic(frame, with_color(frame, orders, colormap); nrow=1) for frame in outs]; fps=5) #hide

# reconstruction by dilation -- `underbuild`
out1 = mreconstruct(dilate, img .- 0.3, img)
# reconstruction by erosion -- `overbuild`
out2 = mreconstruct(erode, img .+ 0.3, img)
mosaic(mask, out1, out2; nrow=1)

mreconstruct(dilate, img .- 0.3, img; dims=1) # reconstruction along colomn
mreconstruct(dilate, img .- 0.3, img, strel_diamond(img)) # diamond shape SE
nothing #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

