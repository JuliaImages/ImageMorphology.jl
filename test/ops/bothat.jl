@testset "bothat" begin
    A = ones(13, 13)
    A[2:3, 2:3] .= 0
    Ae = 1 .- copy(A)
    A[5:9, 5:9] .= 0
    Ao = bothat(A)
    @test Ao == Ae
    @test bothat(A; dims=1:2) == Ao
end
