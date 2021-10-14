@testset "find_boundaries" begin
    A = zeros(Int, 16, 16); A[4:8, 4:8] .= 5; A[4:8, 9:12] .= 6; A[10:12,13:15] .= 3; A[10:12,3:6] .= 9;
    B = A .!= 0
    @test sum(find_boundaries(A)) == 48
    @test hash(find_boundaries(A)) == 0x49da724132e8fcf4
    @test sum(find_boundaries(B)) == 42
    @test hash(find_boundaries(B)) == 0xc6e863a96e69f7ac
end