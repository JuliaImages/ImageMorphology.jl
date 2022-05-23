using ImageMorphology
using ImageCore
using Test
using OffsetArrays
using ImageMetadata
using Suppressor

@test isempty(detect_ambiguities(ImageMorphology))

using Documenter
Base.VERSION >= v"1.6" && doctest(ImageMorphology; manual=false)

include("testutils.jl")
@testset "ImageMorphology" begin
    include("structuring_element.jl")
    include("utils.jl")
    include("extreme_filter.jl")
    include("convexhull.jl")
    include("connected.jl")
    include("ops/dilate.jl")
    include("ops/erode.jl")
    include("ops/opening.jl")
    include("ops/closing.jl")
    include("ops/tophat.jl")
    include("ops/bothat.jl")
    include("ops/mgradient.jl")
    include("ops/morpholaplace.jl")
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
