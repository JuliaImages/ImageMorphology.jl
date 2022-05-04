@testset "label_components" begin
    A = [
        true  true  false true
        true  false true  true
    ]
    lbltarget1 = [
        1 2 0 4
        1 0 3 4
    ]
    @test label_components(A, [1]) == lbltarget1
end

@testset "Erode / dilate" begin
    A = zeros(4, 4, 3)
    A[2, 2, 1] = 0.8
    A[4, 4, 2] = 0.6
    Ad = dilate(A, 1:2)
    Ar = [
        0.8 0.8 0.8 0
        0.8 0.8 0.8 0
        0.8 0.8 0.8 0
        0 0 0 0
    ]
    Ag = [
        0 0 0 0
        0 0 0 0
        0 0 0.6 0.6
        0 0 0.6 0.6
    ]
    @test Ad == cat(Ar, Ag, zeros(4, 4); dims=3)
    @test dilate!(copy(A), 1:2) == Ad
    Ae = erode(Ad, 1:2)
    Ar = [
        0.8 0.8 0 0
        0.8 0.8 0 0
        0 0 0 0
        0 0 0 0
    ]
    Ag = [
        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0.6
    ]
    @test Ae == cat(Ar, Ag, zeros(4, 4); dims=3)
    @test erode!(copy(Ad), 1:2) == Ae
end

@testset "Opening / closing" begin
    A = zeros(4, 4, 3)
    A[2, 2, 1] = 0.8
    A[4, 4, 2] = 0.6
    Ao = opening(A)
    @test opening(A, 1:3) == Ao
    @test opening!(copy(A), 1:3) == Ao
    A = zeros(10, 10)
    A[4:7, 4:7] .= 1
    A[5, 5] = 0
    Ac = closing(A)
    @test closing(A, 1:3) == Ac
    @test closing!(copy(A), 1:3) == Ac
end

@testset "Morphological Top-hat" begin
    A = zeros(13, 13)
    A[2:3, 2:3] .= 1
    A[5:9, 5:9] .= 1
    Ao = tophat(A)
    @test tophat(A, 1:2) == Ao
end

@testset "Morphological Bottom-hat" begin
    A = ones(13, 13)
    A[2:3, 2:3] .= 0
    A[5:9, 5:9] .= 0
    Ao = bothat(A)
    @test bothat(A, 1:2) == Ao
end

@testset "Morphological Gradient" begin
    A = zeros(13, 13)
    A[5:9, 5:9] .= 1
    Ao = morphogradient(A)
    @test morphogradient(A, 1:2) == morphogradient(A)
end

@testset "Morphological Laplacian" begin
    A = zeros(13, 13)
    A[5:9, 5:9] .= 1
    Ao = morpholaplace(A)
    @test morpholaplace(A, 1:2) == Ao
end

@testset "Feature transform/Anisotropic images" begin
    function ind2cart(F)
        s = CartesianIndices(axes(F))
        return map(i -> CartesianIndex(s[i]), F)
    end

    A, w = [false false; false true], (3, 1)
    F = feature_transform(A, w)
    @test F == ind2cart([4 4; 4 4])
    D = distance_transform(F, w)
    @test D ≈ [sqrt(10) 3.0; 1.0 0.0]

    A, w = [false false; false true], (1, 3)
    F = feature_transform(A, w)
    @test F == ind2cart([4 4; 4 4])
    D = distance_transform(F, w)
    @test D ≈ [sqrt(10) 1.0; 3.0 0.0]

    A, w = [true false; false true], (3, 1)
    F = feature_transform(A, w)
    @test F == ind2cart([1 1; 4 4])
    D = distance_transform(F, w)
    @test D == [0 1; 1 0]

    A, w = [true false; false true], (1, 3)
    F = feature_transform(A, w)
    @test F == ind2cart([1 4; 1 4])
    D = distance_transform(F, w)
    @test D == [0 1; 1 0]
end
