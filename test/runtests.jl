using ImageMorphology
using ImageCore
using Test
using OffsetArrays
using OffsetArrays: centered
using ImageMetadata
using Suppressor

@test isempty(detect_ambiguities(ImageMorphology))

using Documenter
Base.VERSION >= v"1.6" && doctest(ImageMorphology, manual = false)

@testset "ImageMorphology" begin
    include("structuring_element.jl")
    include("convexhull.jl")
    include("connected.jl")
    include("dilation_and_erosion.jl")
    include("isboundary.jl")
    include("thinning.jl")
    include("imfill.jl")
    include("maxtree.jl")
    include("feature_transform.jl")
    include("clearborder.jl")
    @info "Beginning deprecation tests, warnings are expected"
    include("deprecations.jl")
end

# multithreaded.jl runs on CI
