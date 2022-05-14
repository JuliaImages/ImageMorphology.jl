@testset "Opening / closing" begin
    A = zeros(4, 4, 3)
    A[2, 2, 1] = 0.8
    A[4, 4, 2] = 0.6
    Ao = opening(A)
    @test Ao == zeros(size(A))
    @test opening(A; dims=1:3) == Ao
    @test opening!(copy(A); dims=1:3) == Ao
    A = zeros(10, 10)
    A[4:7, 4:7] .= 1
    B = copy(A)
    A[5, 5] = 0
    Ac = closing(A)
    @test Ac == B
    @test closing(A; dims=1:3) == Ac
    @test closing!(copy(A); dims=1:3) == Ac
end

@testset "Morphological Top-hat" begin
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

@testset "Morphological Bottom-hat" begin
    A = ones(13, 13)
    A[2:3, 2:3] .= 0
    Ae = 1 .- copy(A)
    A[5:9, 5:9] .= 0
    Ao = bothat(A)
    @test Ao == Ae
    @test bothat(A; dims=1:2) == Ao
end

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

@testset "Morphological Laplacian" begin
    A = zeros(13, 13)
    A[5:9, 5:9] .= 1
    Ao = morpholaplace(A)
    @test morpholaplace(A; dims=1:2) == Ao
    Ae = zeros(13, 13)
    Ae[4:10, 4:10] .= 1
    Ae[5:9, 5:9] .= -1
    Ae[6:8, 6:8] .= 0
    @test Ao == Ae
    Aee = dilate(A) + erode(A) - 2A
    @test Aee == Ae
end
