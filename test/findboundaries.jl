@testset "findboundaries" begin
    # Normal Case 20x20 binary image
    img = [ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        ];
    check_img = Int64.(zeros(size(img)))
    check_img[3, 4:15] .= 1
    check_img[8, 4:15] .= 1
    check_img[3:8,4] .= 1
    check_img[3:8,15] .= 1
    @test check_img == findboundaries(img)

    # background = 1
    check_img = Int64.(zeros(size(img)))
    check_img[2, 3:16] .= 1
    check_img[9, 3:16] .= 1
    check_img[2:9,3] .= 1
    check_img[2:9,16] .= 1
    @test check_img == findboundaries(img, 1)
end