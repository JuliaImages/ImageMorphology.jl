using ImageMorphology
using TestImages
using ImageBase
using ImageShow

img = restrict(Gray.(testimage("blob")))
nothing # hide

out = underbuild(img .- 0.2, img)
mosaic(img, out; nrow=1)

marker, mask = rand(32, 32), rand(32, 32)
complement.(underbuild(marker, mask)) == underbuild(complement.(marker), complement.(mask))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

