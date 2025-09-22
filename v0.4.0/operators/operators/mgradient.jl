using ImageMorphology
using TestImages
using ImageBase
using ImageShow

img = restrict(testimage_dip3e("fig0940")) # rice
img01 = @. Gray(img > 0.5)
mosaic(img, img01; nrow=1)

mosaic(mgradient(img), mgradient(img01); nrow=1)

all(mgradient(rand(32, 32)) .> 0) # always positive

A = rand(32, 32)
mgradient(A) == mgradient(1 .- A) # self-complementary

A = falses(7, 7)
A[3:5, 3:5] .= true
A

Int.(mgradient(A))

Int.(mgradient(A; mode=:external)) # external boundary

Int.(mgradient(A; mode=:internal)) # internal boundary

A = rand(32, 32)
mgradient(A; mode=:external) == mgradient(1 .- A; mode=:internal)

mosaic([mgradient(img01; r) for r in 1:3]; nrow=1)

fout = abs.(ImageBase.FiniteDiff.fdiff(img; dims=1))
mout = mgradient(img; dims=1)
mosaic(fout, mout; nrow=1)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

