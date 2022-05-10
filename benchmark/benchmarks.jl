# Usage:
#     julia benchmark/run_benchmarks.jl
using BenchmarkTools
using ImageCore
using ImageMorphology
using ImageBinarization
using ImageTransformations
using TestImages

on_CI = haskey(ENV, "GITHUB_ACTIONS")

cameraman = Gray{N0f8}.(testimage("cameraman"))
blobs = binarize(Gray.(testimage("blobs")), Otsu()) .> 0.5

tst_sizes = on_CI ? (64, ) : (256, 512)
tst_types = (Gray{N0f8}, Gray{Float32})

const SUITE = BenchmarkGroup()

SUITE["extreme_filter"] = BenchmarkGroup()
let grp = SUITE["extreme_filter"]
    for T in [tst_types..., Gray{Bool}, Int, Bool]
        grp[T] = BenchmarkGroup()
        for sz in tst_sizes
            tst_img = rand(T, sz, sz)

            grp[T]["$sz×$sz"] = BenchmarkGroup()
            se = strel_diamond((3, 3))
            grp[T]["$sz×$sz"]["diamond_r1_fast"] = @benchmarkable extreme_filter(max, $tst_img, $se)
            grp[T]["$sz×$sz"]["diamond_r1_generic"] = @benchmarkable extreme_filter(max, $tst_img, $(collect(se)))
            se = strel_diamond((11, 11))
            grp[T]["$sz×$sz"]["diamond_r5_fast"] = @benchmarkable extreme_filter(max, $tst_img, $se)
            grp[T]["$sz×$sz"]["diamond_r5_generic"] = @benchmarkable extreme_filter(max, $tst_img, $(collect(se)))
        end
    end
end

SUITE["dilatation_and_erosion"] = BenchmarkGroup()
let grp = SUITE["dilatation_and_erosion"]
    grp["erode"] = BenchmarkGroup()
    grp["opening"] = BenchmarkGroup()
    for T in tst_types
        grp["erode"][T] = BenchmarkGroup()
        grp["opening"][T] = BenchmarkGroup()
        for sz in tst_sizes
            tst_img = T.(imresize(T.(cameraman), (sz, sz)))
            grp["erode"][T]["$sz×$sz"] = @benchmarkable erode($tst_img)
            grp["opening"][T]["$sz×$sz"] = @benchmarkable opening($tst_img)
        end
    end
end

SUITE["connected"] = BenchmarkGroup()
let grp = SUITE["connected"]
    grp["label_components"] = @benchmarkable label_components($blobs)
end

SUITE["Maxtree"] = BenchmarkGroup()
let grp = SUITE["Maxtree"]
    grp["area_opening"] = BenchmarkGroup()
    for sz in tst_sizes
        tst_img = (imresize((cameraman), (sz, sz)))
        B = similar(tst_img)
        grp["area_opening"]["$sz×$sz"] = @benchmarkable area_opening($B, min_area=50)
    end
end

SUITE["convexhull"] = BenchmarkGroup()
let grp = SUITE["convexhull"]
    grp["convexhull"] = @benchmarkable convexhull($blobs)
end

SUITE["feature_transform"] = BenchmarkGroup()
let grp = SUITE["feature_transform"]
    grp["feature_transform"] = @benchmarkable feature_transform($blobs)
end
