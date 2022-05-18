@testset "tophat" begin
    A = zeros(13, 13)
    A[2:3, 2:3] .= 1
    Ae = copy(A)
    A[5:9, 5:9] .= 1
    Ao = tophat(A)
    @test Ao == Ae
    @test tophat(A; dims=1:2) == Ae
    Aoo = tophat(Ae)
    @test Aoo == Ae
end
