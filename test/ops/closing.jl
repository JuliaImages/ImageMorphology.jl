@testset "closing" begin
    A = zeros(10, 10)
    A[4:7, 4:7] .= 1
    B = copy(A)
    A[5, 5] = 0
    Ac = closing(A)
    @test Ac == B
    @test closing(A; dims=1:3) == Ac
    @test closing!(similar(A), A, similar(A); dims=1:3) == Ac
end
