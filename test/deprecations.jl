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
