using ImageMorphology
using Test

@testset "ImageMorphology" begin
    include("convexhull.jl")
    include("dilation_and_erosion.jl")
    include("thinning.jl")
end
