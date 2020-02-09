using ImageMorphology
using Test
using OffsetArrays

@testset "ImageMorphology" begin
    include("convexhull.jl")
    include("connected.jl")
    include("dilation_and_erosion.jl")
    include("thinning.jl")
    include("imfill.jl")
    include("maxtree.jl")
end
