opening_ref(img, dims, r) = opening_ref(img, strel_box(img, dims; r))
opening_ref(img, se) = dilate(erode(img, se), se)

@testset "opening" begin
    test_types = [Bool, Int, Float64, Gray{N0f8}, Gray{Float64}]
    for T in test_types
        for N in (1, 2, 3)
            sz = ntuple(_ -> 32, N)
            img = rand(T, sz...)

            out = opening(img)
            @test out == opening(img, strel_box(img, ntuple(identity, N); r=1))
            @test out == opening_ref(img, ntuple(identity, N), 1)

            @test opening(img; dims=(1,), r=2) == opening_ref(img, (1,), 2)

            se = centered(rand(Bool, ntuple(_ -> 3, N)))
            @test opening(img, se) == opening_ref(img, se)
        end
    end

    img = rand(1:5, 7, 7)
    out = similar(img)
    opening!(out, img, similar(img))
    @test out == opening(img)

    img = rand(RGB, 7, 7)
    @test_throws ArgumentError("color image is not supported") opening(img)
end
