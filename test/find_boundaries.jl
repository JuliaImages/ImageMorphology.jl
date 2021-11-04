import ImageMorphology: find_boundaries_thick, find_boundaries_dilate_erode

@testset "find_boundaries" begin
    A = zeros(Int, 16, 16); A[4:8, 4:8] .= 5; A[4:8, 9:12] .= 6; A[10:12,13:15] .= 3; A[10:12,3:6] .= 9;
    OA = OffsetArray(A, -1, -1)
    FA = Float32.(A./2)
    B = A .!= 0
    OB = OffsetArray(B, -1, -1)
    GB = Gray{Bool}.(B)
    normA = A ./ maximum(A)
    C = copy(A)
    C[A .== 0] .= 1
    normC = C ./ maximum(C)
    FC = Float32.(C./2)
    @test sum(find_boundaries(A)) == 48
    @test sum(find_boundaries(OA)) == 48
    @test sum(find_boundaries(FA)) === 48.0f0
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
    @test sum(find_boundaries(FC)) === 107.0f0
    @test sum(find_boundaries(FC, background = 1/2)) === 48.0f0
    @test sum(find_boundaries(Gray{Float32}.(normC), background = 0)) == 107
    @test sum(find_boundaries(Gray{Float32}.(normC), background = Gray{Float32}(1/9))) == 48
    @test sum(find_boundaries(RGB{N0f8}.(normC), background = RGB{N0f8}(1/9))) .== RGB{Float64}(48, 48, 48)

    # Tests from pull request #32
    # Normal Case 20x20 binary image
    img = [ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        ];
    check_img = Int64.(zeros(size(img)))
    check_img[3, 4:15] .= 1
    check_img[8, 4:15] .= 1
    check_img[3:8,4] .= 1
    check_img[3:8,15] .= 1
    @test check_img == find_boundaries(img)

    # background = 1
    check_img = Int64.(zeros(size(img)))
    check_img[2, 3:16] .= 1
    check_img[9, 3:16] .= 1
    check_img[2:9,3] .= 1
    check_img[2:9,16] .= 1
    @test check_img == find_boundaries(img; background = 1)
end

@testset "find_boundaries_thick" begin
    A = zeros(Int, 16, 16); A[4:8, 4:8] .= 5; A[4:8, 9:12] .= 6; A[10:12,13:15] .= 3; A[10:12,3:6] .= 9;
    OA = OffsetArray(A, -1, -1)
    FA = Float32.(A./2)
    B = A .!= 0
    OB = OffsetArray(B, -1, -1)
    GB = Gray{Bool}.(B)
    normA = A ./ maximum(A)
    C = copy(A)
    C[A .== 0] .= 1
    FC = Float32.(C./2)
    @test sum(find_boundaries_thick(A)) == 107
    @test sum(find_boundaries_thick(OA)) == 107
    @test sum(find_boundaries_thick(FA)) == 107
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
    @test sum(find_boundaries_thick(C)) == 107
    @test sum(find_boundaries_thick(FC)) == 107
    img = [ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        ];
    check_img = [
         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
        0  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  0  0  0
        0  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  0  0  0
        0  0  1  1  0  0  0  0  0  0  0  0  0  0  1  1  0  0  0  0
        0  0  1  1  0  0  0  0  0  0  0  0  0  0  1  1  0  0  0  0
        0  0  1  1  0  0  0  0  0  0  0  0  0  0  1  1  0  0  0  0
        0  0  1  1  0  0  0  0  0  0  0  0  0  0  1  1  0  0  0  0
        0  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  0  0  0
        0  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  0  0  0
        0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
    ]
    @test find_boundaries_thick(img) == check_img
end

@testset "find_boundaries == find_boundaries_thick .& image" begin
    A = zeros(Int, 16, 16); A[4:8, 4:8] .= 5; A[4:8, 9:12] .= 6; A[10:12,13:15] .= 3; A[10:12,3:6] .= 9;
    B = A .!= 0
    C = copy(A)
    C[A .== 0] .= 1
    @test find_boundaries(A) == find_boundaries_thick(A) .& (A .!= 0)
    @test find_boundaries(B) == find_boundaries_thick(B) .& B
    @test find_boundaries(B, background=false) == find_boundaries_thick(B) .& B
    @test find_boundaries(B, background=true) == find_boundaries_thick(B) .& .!B
    @test find_boundaries(C) == find_boundaries_thick(C)
    @test find_boundaries(C, background = 1) == find_boundaries_thick(C) .& (C .!= 1)
    @test find_boundaries(A) == find_boundaries_dilate_erode(A)
end