@testset "clearborder" begin
    #Case when given border width is more than image size
    img = [1 0 1 1 1 1
            0 1 1 1 0 0
            1 1 0 0 0 1
            0 1 0 1 0 1
            1 1 0 0 0 1
            0 0 1 1 0 0]
    @test_throws ArgumentError clearborder(img,7)

    #Normal Case
    img = [0 0 0 0 0 0 0 1 0
            0 0 0 0 1 0 0 0 0
            1 0 0 1 0 1 0 0 0
            0 0 1 1 1 1 1 0 0
            0 1 1 1 1 1 1 1 0
            0 0 0 0 0 0 0 0 0]
    cleared_img = clearborder(img)
    check_img = copy(img)
    check_img[3,1] = 0
    check_img[1,8] = 0
    @test cleared_img == check_img

    cleared_img = clearborder(img,2)
    @test cleared_img == fill!(similar(img), zero(eltype(img)))

    cleared_img = clearborder(img,2,10)
    @test cleared_img == 10*fill!(similar(img), one(eltype(img)))

    #Multidimentional Case
    img = cat([0 0 0 0;
                    0 0 0 0;
                    0 0 0 0;
                    1 0 0 0],
                [0 0 0 0;
                    0 1 1 0;
                    0 0 1 0;
                    0 0 0 0],
                [0 0 0 0;
                    0 0 0 0;
                    0 0 0 0;
                    0 0 0 0], dims=3)
    cleared_img = clearborder(img)
    check_img = copy(img)
    check_img[4,1,1] = 0
    @test cleared_img == check_img

    cleared_img = clearborder(img,2)
    @test cleared_img == fill!(similar(img), zero(eltype(img)))

    cleared_img = clearborder(img,2,10)
    @test cleared_img == 10*fill!(similar(img), one(eltype(img)))

    #Grayscale input image Case
    img = [1 2 3 1 2
            3 3 5 4 2
            3 4 5 4 2
            3 3 2 1 2]
    cleared_img = clearborder(img)
    check_img = [0 0 0 0 0
                    0 0 5 4 0
                    0 4 5 4 0
                    0 0 0 0 0]
    @test cleared_img == check_img
end