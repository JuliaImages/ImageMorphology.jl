@testset "opening" begin
    A = zeros(4, 4, 3)
    A[2, 2, 1] = 0.8
    A[4, 4, 2] = 0.6
    Ao = opening(A)
    @test Ao == zeros(size(A))
    @test opening(A; dims=1:3) == Ao
    @test opening!(similar(A), A, similar(A); dims=1:3) == Ao
end
