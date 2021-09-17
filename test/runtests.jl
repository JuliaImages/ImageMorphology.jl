using ImageMorphology
using ImageCore
using Test
using OffsetArrays
using ImageMetadata

@test isempty(detect_ambiguities(ImageMorphology))

using Documenter
Base.VERSION >= v"1.6" && doctest(ImageMorphology, manual = false)

@testset "ImageMorphology" begin
    include("convexhull.jl")
    include("connected.jl")
    include("dilation_and_erosion.jl")
    include("thinning.jl")
    include("imfill.jl")
    include("maxtree.jl")
    include("bwdist.jl")
    include("clearborder.jl")
end
