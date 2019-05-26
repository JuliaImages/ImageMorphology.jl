using ImageMorphology
using Images
using Test

@testset "ImageMorphology" begin
    include("convexhull.jl")
    include("connected.jl")
    include("dilation_and_erosion.jl")
    include("thinning.jl")
    include("imfill.jl")
end
