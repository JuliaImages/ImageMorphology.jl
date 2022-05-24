beucher_gradient_ref(img, dims, r) = beucher_gradient_ref(img, strel_box(img, dims; r))
beucher_gradient_ref(img, se) = float.(dilate(img, se)) .- float.(erode(img, se))
external_gradient_ref(img, dims, r) = external_gradient_ref(img, strel_box(img, dims; r))
external_gradient_ref(img, se) = float.(dilate(img, se)) .- img
internal_gradient_ref(img, dims, r) = internal_gradient_ref(img, strel_box(img, dims; r))
internal_gradient_ref(img, se) = img .- float.(erode(img, se))

@testset "mgradient" begin
    test_ranges = [Bool, 1:10, Float64, Gray{N0f8}, Gray{Float64}]
    modes = [
        (:beucher, beucher_gradient_ref),
        (:external, external_gradient_ref),
        (:internal, internal_gradient_ref),
    ]
    for (mode, ref_fun) in modes
        for T in test_ranges
            for N in (1, 2, 3)
                sz = ntuple(_ -> 32, N)
                img = rand(T, sz...)

                out = mgradient(img; mode)
                @test out == mgradient(img, strel_box(img, ntuple(identity, N); r=1); mode)
                @test out == ref_fun(img, ntuple(identity, N), 1)

                @test mgradient(img; dims=(1,), r=2, mode) == ref_fun(img, (1,), 2)

                se = rand_se_mask(3, N; symmetric=true)
                @test mgradient(img, se; mode) == ref_fun(img, se)
            end
        end
    end

    img = rand(1:5, 7, 7)
    out = similar(img)
    mgradient!(out, img, similar(img))
    @test out == mgradient(img)

    img = rand(RGB, 7, 7)
    @test_throws ArgumentError("color image is not supported") mgradient(img)

    img = rand(1:10, 7, 7)
    se = Bool[1 0 0; 0 1 0; 0 0 0]
    msg = "structuring element must be symmetric with respect to its center"
    @test_throws ArgumentError(msg) mlaplacian(img, se)
end
