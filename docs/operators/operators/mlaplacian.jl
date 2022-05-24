# ---
# title: "`mlaplacian`"
# id: op_mlaplacian
# cover: assets/mlaplacian.png
# ---

# Laplacian operator is defined as the difference between external gradient and internal
# gradient -- `mgradient(img; mode=:external) - mgradient(img; mode=:internal)`.

using ImageMorphology
using TestImages
using ImageBase
using ImageShow

img = restrict(testimage_dip3e("fig0940")) # rice
img01 = @. Gray(img > 0.5)
mosaic(img, img01; nrow=1)

#-

mosaic(mlaplacian(img), mlaplacian(img01); nrow=1)

# Note that laplacian can produce negative values:

A = falses(7, 7)
A[3:5, 3:5] .= true
A[4, 4] = false
Int.(mlaplacian(A))

# ## See also

# ImageBase.jl provides the finite-difference version of laplace operator `flaplacian`.
# ImageFiltering.jl provides `Laplacian` and `LaplacianOfGaussian` kernels.

fout = abs.(ImageBase.FiniteDiff.flaplacian(img))
mout = abs.(mlaplacian(img))
mosaic(fout, mout; nrow=1)

# For a comprehensive and more accurate documentation, please check the [`mlaplacian`](@ref)
# reference page.

# save cover image #src
using FileIO #src
mkpath("assets") #src
img = Gray.(restrict(testimage("blobs"))) #src
save("assets/mlaplacian.png", abs.(mlaplacian(img; r=2))) #src
