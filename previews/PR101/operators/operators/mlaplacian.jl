using ImageMorphology
using TestImages
using ImageBase
using ImageShow

img = restrict(testimage_dip3e("fig0940")) # rice
img01 = @. Gray(img > 0.5)
mosaic(img, img01; nrow=1)

mosaic(mlaplacian(img), mlaplacian(img01); nrow=1)

A = falses(7, 7)
A[3:5, 3:5] .= true
A[4, 4] = false
Int.(mlaplacian(A))

fout = abs.(ImageBase.FiniteDiff.flaplacian(img))
mout = abs.(mlaplacian(img))
mosaic(fout, mout; nrow=1)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

