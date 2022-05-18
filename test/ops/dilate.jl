@testset "dilate" begin
    A = zeros(4, 4, 3)
    A[2, 2, 1] = 0.8
    A[4, 4, 2] = 0.6
    Ae = erode(A)
    @test Ae == zeros(size(A))
    Ad = dilate(A; dims=1:2)
    Ar = [
        0.8 0.8 0.8 0
        0.8 0.8 0.8 0
        0.8 0.8 0.8 0
        0 0 0 0
    ]
    Ag = [
        0 0 0 0
        0 0 0 0
        0 0 0.6 0.6
        0 0 0.6 0.6
    ]
    @test Ad == cat(Ar, Ag, zeros(4, 4); dims=3)
    @test dilate!(similar(A), A; dims=1:2) == Ad
    Ae = erode(Ad; dims=1:2)
    Ar = [
        0.8 0.8 0 0
        0.8 0.8 0 0
        0 0 0 0
        0 0 0 0
    ]
    Ag = [
        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0.6
    ]
    @test Ae == cat(Ar, Ag, zeros(4, 4); dims=3)
    @test erode!(similar(Ad), Ad; dims=1:2) == Ae
    # issue Images.jl #311
    @test dilate(trues(3)) == trues(3)
    # ImageMeta
    @test arraydata(dilate(ImageMeta(A))) == dilate(A)
    @test arraydata(dilate(ImageMeta(A); dims=1:2)) == dilate(A; dims=1:2)
    @test arraydata(erode(ImageMeta(A))) == erode(A)
    @test arraydata(erode(ImageMeta(A); dims=1:2)) == erode(A; dims=1:2)
end
