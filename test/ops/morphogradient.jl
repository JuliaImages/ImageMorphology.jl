beucher_gradient_ref(img, dims, r) = beucher_gradient_ref(img, strel_box(img, dims; r))
beucher_gradient_ref(img, se) = float.(dilate(img, se)) .- float.(erode(img, se))

@testset "morphogradient" begin
    test_ranges = [Bool, 1:10, Float64, Gray{N0f8}, Gray{Float64}]
    for T in test_ranges
        for N in (1, 2, 3)
            sz = ntuple(_ -> 32, N)
            img = rand(T, sz...)

            out = morphogradient(img)
            @test out == morphogradient(img, strel_box(img, ntuple(identity, N); r=1))
            @test out == beucher_gradient_ref(img, ntuple(identity, N), 1)

            @test morphogradient(img; dims=(1,), r=2) == beucher_gradient_ref(img, (1,), 2)

            se = rand_se_mask(3, N; symmetric=true)
            @test morphogradient(img, se) == beucher_gradient_ref(img, se)
        end
    end

    # TODO(johnnychen94): support this
    # img = rand(1:5, 7, 7)
    # out = similar(img)
    # morphogradient!(out, img, similar(img))
    # @test out == morphogradient(img)

    img = rand(RGB, 7, 7)
    @test_throws ArgumentError("color image is not supported") morphogradient(img)

    img = rand(1:10, 7, 7)
    se = Bool[1 0 0; 0 1 0; 0 0 0]
    msg = "structuring element must be symmetric with respect to its center"
    @test_throws ArgumentError(msg) morpholaplace(img, se)
end
