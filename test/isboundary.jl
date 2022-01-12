import ImageMorphology: isboundary_thick, isboundary_dilate_erode

@testset "isboundary" begin
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
    @test sum(isboundary(A)) == 48
    @test sum(isboundary(OA)) == 48
    @test sum(isboundary(FA)) === 48.0f0
    @test sum(isboundary(Gray{Float32}.(normA))) == 48
    @test sum(isboundary(RGB{N0f8}.(normA))) .== RGB{Float64}(48, 48, 48)
    @test sum(isboundary(A; dims = 1)) == 32
    @test sum(isboundary(A; dims = 2)) == 32
    @test sum(isboundary(GB)) == 42
    @test sum(isboundary(B)) == 42
    @test sum(isboundary(OB)) == 42
    @test sum(isboundary(B; dims = 1)) == 32
    @test sum(isboundary(B; dims = 2)) == 22
    @test sum(isboundary(C)) == 107
    @test sum(isboundary(C, background = 1)) == 48
    @test sum(isboundary(FC)) === 107.0f0
    @test sum(isboundary(FC, background = 1/2)) === 48.0f0
    @test sum(isboundary(Gray{Float32}.(normC), background = 0)) == 107
    @test sum(isboundary(Gray{Float32}.(normC), background = Gray{Float32}(1/9))) == 48
    @test sum(isboundary(RGB{N0f8}.(normC), background = RGB{N0f8}(1/9))) .== RGB{Float64}(48, 48, 48)

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
    @test check_img == isboundary(img)

    # background = 1
    check_img = Int64.(zeros(size(img)))
    check_img[2, 3:16] .= 1
    check_img[9, 3:16] .= 1
    check_img[2:9,3] .= 1
    check_img[2:9,16] .= 1
    @test check_img == isboundary(img; background = 1)
end

@testset "isboundary_thick" begin
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
    @test sum(isboundary_thick(A)) == 107
    @test sum(isboundary_thick(OA)) == 107
    @test sum(isboundary_thick(FA)) == 107
    @test sum(isboundary_thick(Gray{Float32}.(normA))) == 107
    # MethodError: no method matching isless(::RGB{N0f8}, ::RGB{N0f8})
    @test_broken sum(isboundary_broken(RGB{N0f8}.(normA))) .== RGB{Float64}(107, 107, 107)
    @test sum(isboundary_thick(A; dims = 1)) == 61
    @test sum(isboundary_thick(A; dims = 2)) == 54
    @test sum(isboundary_thick(GB)) == 101
    @test sum(isboundary_thick(B)) == 101
    @test sum(isboundary_thick(OB)) == 101
    @test sum(isboundary_thick(B; dims = 1)) == 61
    @test sum(isboundary_thick(B; dims = 2)) == 44
    @test sum(isboundary_thick(C)) == 107
    @test sum(isboundary_thick(FC)) == 107
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
    @test isboundary_thick(img) == check_img
end

@testset "isboundary == isboundary_thick .& image" begin
    A = zeros(Int, 16, 16); A[4:8, 4:8] .= 5; A[4:8, 9:12] .= 6; A[10:12,13:15] .= 3; A[10:12,3:6] .= 9;
    B = A .!= 0
    C = copy(A)
    C[A .== 0] .= 1
    @test isboundary(A) == isboundary_thick(A) .& (A .!= 0)
    @test isboundary(B) == isboundary_thick(B) .& B
    @test isboundary(B, background=false) == isboundary_thick(B) .& B
    @test isboundary(B, background=true) == isboundary_thick(B) .& .!B
    @test isboundary(C) == isboundary_thick(C)
    @test isboundary(C, background = 1) == isboundary_thick(C) .& (C .!= 1)
    @test isboundary(A) == isboundary_dilate_erode(A)
end