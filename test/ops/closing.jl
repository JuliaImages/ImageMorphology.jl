closing_ref(img, dims, r) = closing_ref(img, strel_box(img, dims; r))
closing_ref(img, se) = erode(dilate(img, se), se)

@testset "closing" begin
    test_types = [Bool, Int, Float64, Gray{N0f8}, Gray{Float64}]
    for T in test_types
        for N in (1, 2, 3)
            sz = ntuple(_ -> 32, N)
            img = rand(T, sz...)

            out = closing(img)
            @test out == closing(img, strel_box(img, ntuple(identity, N); r=1))
            @test out == closing_ref(img, ntuple(identity, N), 1)

            @test closing(img; dims=(1,), r=2) == closing_ref(img, (1,), 2)

            se = centered(rand(Bool, ntuple(_ -> 3, N)))
            @test closing(img, se) == closing_ref(img, se)
        end
    end

    img = rand(1:5, 7, 7)
    out = similar(img)
    closing!(out, img, similar(img))
    @test out == closing(img)

    img = rand(RGB, 7, 7)
    @test_throws ArgumentError("color image is not supported") closing(img)
end
