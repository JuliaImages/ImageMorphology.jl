@testset "Thinning" begin
    # for □ like figure
    img = falses(20, 20)
    for i in 7:(7 + 4), j in 7:(7 + 4)
        img[i, j] = img[j, i] = true
    end
    thin = thinning(img)
    @test count(thin .== 1) == 1
    @test thin[9, 9] == 1

    # for + like figure
    img = falses(20, 20)
    for i in 8:(8 + 5), j in 4:(13 + 4)
        img[i, j] = true
        img[j, i] = true
    end
    ans = falses(size(img))
    ans[[111, 131, 151, 171, 191, 192, 193, 194, 195, 211, 231, 251, 271]] .= 1
    @test thinning(img) == ans

    img = Bool[
        1 1 1 1 1 1
        0 0 0 0 0 0
        1 1 1 1 1 1
        1 1 1 1 1 1
        0 0 1 1 0 0
        0 0 0 0 0 0
    ]
    ans = Bool[
        1 1 1 1 1 1
        0 0 0 0 0 0
        0 0 0 0 0 0
        1 1 1 1 1 0
        0 0 0 0 0 0
        0 0 0 0 0 0
    ]
    @test thinning(img) == ans

    # already thinned
    img = Bool[
        0 0 0 0 0 0
        0 0 0 1 1 0
        0 0 1 0 0 0
        0 0 1 0 0 0
        0 0 0 1 0 0
        0 0 0 0 0 0
        0 0 0 0 0 0
    ]
    @test thinning(img) == img
end
