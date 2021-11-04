import ImageMorphology: find_boundaries_thick

@testset "find_boundaries" begin
    A = zeros(Int, 16, 16); A[4:8, 4:8] .= 5; A[4:8, 9:12] .= 6; A[10:12,13:15] .= 3; A[10:12,3:6] .= 9;
    OA = OffsetArray(A, -1, -1)
    B = A .!= 0
    OB = OffsetArray(B, -1, -1)
    GB = Gray{Bool}.(B)
    normA = A ./ maximum(A)
    C = copy(A)
    C[A .== 0] .= 1
    normC = C ./ maximum(C)
    @test sum(find_boundaries(A)) == 48
    @test sum(find_boundaries(OA)) == 48
    @test sum(find_boundaries(Gray{Float32}.(normA))) == 48
    @test sum(find_boundaries(RGB{N0f8}.(normA))) .== RGB{Float64}(48, 48, 48)
    @test sum(find_boundaries(A; dims = 1)) == 32
    @test sum(find_boundaries(A; dims = 2)) == 32
    @test sum(find_boundaries(GB)) == 42
    @test sum(find_boundaries(B)) == 42
    @test sum(find_boundaries(OB)) == 42
    @test sum(find_boundaries(B; dims = 1)) == 32
    @test sum(find_boundaries(B; dims = 2)) == 22
    @test sum(find_boundaries(C)) == 107
    @test sum(find_boundaries(C, background = 1)) == 48
    @test sum(find_boundaries(Gray{Float32}.(normC), background = 0)) == 107
    @test sum(find_boundaries(Gray{Float32}.(normC), background = Gray{Float32}(1/9))) == 48
    @test sum(find_boundaries(RGB{N0f8}.(normC), background = RGB{N0f8}(1/9))) .== RGB{Float64}(48, 48, 48)
end

@testset "find_boundaries_thick" begin
    A = zeros(Int, 16, 16); A[4:8, 4:8] .= 5; A[4:8, 9:12] .= 6; A[10:12,13:15] .= 3; A[10:12,3:6] .= 9;
    OA = OffsetArray(A, -1, -1)
    B = A .!= 0
    OB = OffsetArray(B, -1, -1)
    GB = Gray{Bool}.(B)
    normA = A ./ maximum(A)
    C = copy(A)
    C[A .== 0] .= 1
    @test sum(find_boundaries_thick(A)) == 107
    @test sum(find_boundaries_thick(OA)) == 107
    @test sum(find_boundaries_thick(Gray{Float32}.(normA))) == 107
    # MethodError: no method matching isless(::RGB{N0f8}, ::RGB{N0f8})
    @test_broken sum(find_boundaries_broken(RGB{N0f8}.(normA))) .== RGB{Float64}(107, 107, 107)
    @test sum(find_boundaries_thick(A; dims = 1)) == 61
    @test sum(find_boundaries_thick(A; dims = 2)) == 54
    @test sum(find_boundaries_thick(GB)) == 101
    @test sum(find_boundaries_thick(B)) == 101
    @test sum(find_boundaries_thick(OB)) == 101
    @test sum(find_boundaries_thick(B; dims = 1)) == 61
    @test sum(find_boundaries_thick(B; dims = 2)) == 44
end

@testset "find_boundaries == find_boundaries_thick .& image" begin
    A = zeros(Int, 16, 16); A[4:8, 4:8] .= 5; A[4:8, 9:12] .= 6; A[10:12,13:15] .= 3; A[10:12,3:6] .= 9;
    B = A .!= 0
    @test find_boundaries(A) == find_boundaries_thick(A) .& (A .!= 0)
    @test find_boundaries(B) == find_boundaries_thick(B) .& B
    @test find_boundaries(B, background=false) == find_boundaries_thick(B) .& B
    @test find_boundaries(B, background=true) == find_boundaries_thick(B) .& .!B
end