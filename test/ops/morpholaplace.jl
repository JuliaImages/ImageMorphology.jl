morpholaplace_ref(img, dims, r) = morpholaplace_ref(img, strel_box(img, dims; r))
morpholaplace_ref(img, se) = float.(dilate(img, se)) .+ float.(erode(img, se)) .- 2 .* float.(img)

@testset "morpholaplace" begin
    test_ranges = [Bool, 1:10, Float64, Gray{N0f8}, Gray{Float64}]
    for T in test_ranges
        for N in (1, 2, 3)
            sz = ntuple(_ -> 32, N)
            img = rand(T, sz...)

            out = morpholaplace(img)
            @test out == morpholaplace(img, strel_box(img, ntuple(identity, N); r=1))
            @test out == morpholaplace_ref(img, ntuple(identity, N), 1)

            @test morpholaplace(img; dims=(1,), r=2) == morpholaplace_ref(img, (1,), 2)

            se = centered(rand(Bool, ntuple(_ -> 3, N)))
            @test morpholaplace(img, se) == morpholaplace_ref(img, se)
        end
    end

    # TODO(johnnychen94): support this
    # img = rand(1:5, 7, 7)
    # out = similar(img)
    # morpholaplace!(out, img, similar(img))
    # @test out == morpholaplace(img)

    img = rand(RGB, 7, 7)
    @test_throws ArgumentError("color image is not supported") morpholaplace(img)
end
