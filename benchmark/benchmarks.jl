Usage:
    julia benchmark/run_benchmarks.jl
using BenchmarkTools
using ImageCore
using ImageMorphology
using ImageBinarization
using ImageTransformations
using ImageFiltering
using TestImages

on_CI = haskey(ENV, "GITHUB_ACTIONS")

cameraman = Gray{N0f8}.(testimage("cameraman"))
blobs = binarize(Gray.(testimage("blobs")), Otsu()) .> 0.5

tst_sizes = on_CI ? (256,) : (256, 512)
tst_types = (Gray{N0f8}, Gray{Float32})

const SUITE = BenchmarkGroup()

SUITE["extreme_filter"] = BenchmarkGroup()
let grp = SUITE["extreme_filter"]
    for T in [tst_types..., Int]
        grp[T] = BenchmarkGroup()
        for sz in tst_sizes
            tst_img = rand(T, sz, sz)

            grp[T]["$sz×$sz"] = BenchmarkGroup()
            se = strel_diamond((3, 3))
            grp[T]["$sz×$sz"]["r1_diamond"] = @benchmarkable extreme_filter(max, $tst_img, $se)
            se = strel_box((3, 3))
            grp[T]["$sz×$sz"]["r1_box"] = @benchmarkable extreme_filter(max, $tst_img, $se)
            se = centered(collect(se))
            grp[T]["$sz×$sz"]["r1_generic"] = @benchmarkable extreme_filter(max, $tst_img, $se)
            se = strel_diamond((11, 11))
            grp[T]["$sz×$sz"]["r5_diamond"] = @benchmarkable extreme_filter(max, $tst_img, $se)
            se = strel_box((11, 11))
            grp[T]["$sz×$sz"]["r5_box"] = @benchmarkable extreme_filter(max, $tst_img, $se)
            se = centered(collect(se))
            grp[T]["$sz×$sz"]["r5_generic"] = @benchmarkable extreme_filter(max, $tst_img, $se)
        end
    end
end
let grp = SUITE["extreme_filter"]
    T = Bool
    grp[T] = BenchmarkGroup()
    for sz in tst_sizes
        grp[T]["$sz×$sz"] = BenchmarkGroup()
        for (cname, cr) in [("worst", 0.95), ("best", 0.05), ("random", 0.5)]
            tst_img = fill(zero(T), sz, sz)
            tst_img[rand(sz, sz) .>= cr] .= true

            se = strel_diamond((3, 3))
            grp[T]["$sz×$sz"]["r1_diamond_$cname"] = @benchmarkable extreme_filter(max, $tst_img, $se)
            se = centered(collect(se))
            grp[T]["$sz×$sz"]["r1_bool_$cname"] = @benchmarkable extreme_filter(max, $tst_img, $se)

            se = strel_diamond((11, 11))
            grp[T]["$sz×$sz"]["r5_diamond_$cname"] = @benchmarkable extreme_filter(max, $tst_img, $se)
            se = centered(collect(se))
            grp[T]["$sz×$sz"]["r5_bool_$cname"] = @benchmarkable extreme_filter(max, $tst_img, $se)
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

SUITE["mreconstruct"] = BenchmarkGroup()
let grp = SUITE["mreconstruct"]
    @assert eltype(blobs) == Bool
    erode_blobs = erode(blobs)
    grp["Bool"] = @benchmarkable mreconstruct(dilate, $erode_blobs, $blobs)
    for T in tst_types
        grp["$T"] = @benchmarkable mreconstruct(dilate, $(T.(erode_blobs)), $(T.(blobs)))
    end
end

SUITE["connected"] = BenchmarkGroup()
let grp = SUITE["connected"]
    grp["label_components"] = @benchmarkable label_components($blobs)
    grp["label_flatzones"] = @benchmarkable label_flatzones($cameraman, trues(3, 3))
    grp["label_lambdaflatzones"] = @benchmarkable label_lambdaflatzones($cameraman, trues(3, 3), Gray{N0f8}(1.0 / 255.0))
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

SUITE["leveling"] = BenchmarkGroup()
let grp = SUITE["leveling"]
    gaus = Gray{N0f8}.(imfilter(cameraman, Kernel.gaussian(5)))
    grp["low_leveling"] = @benchmarkable low_leveling(cameraman, $gaus)
    grp["high_leveling"] = @benchmarkable high_leveling(cameraman, $gaus)
    grp["leveling"] = @benchmarkable leveling(cameraman, $gaus)
end

SUITE["extremum"] = BenchmarkGroup()
let grp = SUITE["extremum"]
    grp["hmaxima"] = BenchmarkGroup()
    for sz in tst_sizes
        tst_img = reinterpret(N0f8, imresize((cameraman), (sz, sz)))
        grp["hmaxima"]["$sz×$sz"] = @benchmarkable hmaxima($tst_img, reinterpret(N0f8, (UInt8)(3)))
    end
    grp["regional_maxima"] = BenchmarkGroup()
    for sz in tst_sizes
        tst_img = reinterpret(N0f8, imresize((cameraman), (sz, sz)))
        grp["regional_maxima"]["$sz×$sz"] = @benchmarkable regional_maxima($tst_img)
    end
end

SUITE["clearborder"] = BenchmarkGroup()
let grp = SUITE["clearborder"]
    grp["clearborder"] = @benchmarkable clearborder($blobs)
end

SUITE["fillhole"] = BenchmarkGroup()
let grp = SUITE["fillhole"]
    grp["fillhole"] = @benchmarkable fillhole($blobs)
end