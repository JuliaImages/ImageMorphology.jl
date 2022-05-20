tophat_ref(img, dims, r) = tophat_ref(img, strel_box(img, dims; r))
tophat_ref(img, se) = float.(img) .- float.(opening(img, se))

@testset "tophat" begin
    test_ranges = [Bool, 1:10, Float64, Gray{N0f8}, Gray{Float64}]
    for T in test_ranges
        for N in (1, 2, 3)
            sz = ntuple(_ -> 32, N)
            img = rand(T, sz...)

            out = tophat(img)
            @test out == tophat(img, strel_box(img, ntuple(identity, N); r=1))
            @test out == tophat_ref(img, ntuple(identity, N), 1)

            @test tophat(img; dims=(1,), r=2) == tophat_ref(img, (1,), 2)

            se = centered(rand(Bool, ntuple(_ -> 3, N)))
            @test tophat(img, se) == tophat_ref(img, se)
        end
    end

    img = rand(1:5, 7, 7)
    out = similar(img)
    tophat!(out, img, similar(img))
    @test out == tophat(img)

    img = rand(RGB, 7, 7)
    @test_throws ArgumentError("color image is not supported") tophat(img)
end
