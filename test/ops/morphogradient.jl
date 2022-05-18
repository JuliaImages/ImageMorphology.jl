@testset "Morphological Gradient" begin
    A = zeros(13, 13)
    A[5:9, 5:9] .= 1
    Ao = morphogradient(A)
    @test morphogradient(A; dims=1:2) == morphogradient(A)
    Ae = zeros(13, 13)
    Ae[4:10, 4:10] .= 1
    Ae[6:8, 6:8] .= 0
    @test Ao == Ae
    Aee = dilate(A) - erode(A)
    @test Aee == Ae
end