using ImageMorphology
using ImageCore
using ImageCore.PaddedViews
using Test
using OffsetArrays
using ImageMetadata
using Suppressor

@test isempty(detect_ambiguities(ImageMorphology))

if Base.VERSION >= v"1.8"
	using Documenter
    doctest(ImageMorphology; manual=false)
end

include("testutils.jl")

include("structuring_element.jl")
include("utils.jl")
include("extreme_filter.jl")
include("extremum.jl")
include("convexhull.jl")
include("connected.jl")
include("ops/dilate.jl")
include("ops/erode.jl")
include("ops/opening.jl")
include("ops/closing.jl")
include("ops/tophat.jl")
include("ops/bothat.jl")
include("ops/mgradient.jl")
include("ops/mlaplacian.jl")
include("ops/mreconstruct.jl")
include("ops/underbuild.jl")
include("ops/overbuild.jl")
include("isboundary.jl")
include("thinning.jl")
include("imfill.jl")
include("maxtree.jl")
include("feature_transform.jl")
include("leveling.jl")
include("clearborder.jl")

@info "Beginning deprecation tests, warnings are expected"
include("deprecations.jl")

# multithreaded.jl runs on CI
