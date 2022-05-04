@testset "Label components" begin
    A = [
        true  true  false true
        true  false true  true
    ]
    lbltarget = [
        1 1 0 2
        1 0 2 2
    ]
    lbltarget1 = [
        1 2 0 4
        1 0 3 4
    ]
    @test label_components(A) == lbltarget
    @test label_components(A; dims=:) == lbltarget
    @test label_components(A; dims=1) == lbltarget1
    connectivity = [
        false true  false
        true  false true
        false true  false
    ]
    @test label_components(A, connectivity) == lbltarget
    connectivity = trues(3, 3)
    lbltarget2 = [
        1 1 0 1
        1 0 1 1
    ]
    @test label_components(A, connectivity) == lbltarget2
    @test component_boxes(lbltarget) ==
        Vector{Tuple}[[(1, 2), (2, 3)], [(1, 1), (2, 2)], [(1, 3), (2, 4)]]
    @test component_lengths(lbltarget) == [2, 3, 3]
    @test component_indices(lbltarget) == Array{Int}[[4, 5], [1, 2, 3], [6, 7, 8]]
    @test component_subscripts(lbltarget) ==
        Array{Tuple}[[(2, 2), (1, 3)], [(1, 1), (2, 1), (1, 2)], [(2, 3), (1, 4), (2, 4)]]
    @test @inferred(component_centroids(lbltarget)) ==
        Tuple[(1.5, 2.5), (4 / 3, 4 / 3), (5 / 3, 11 / 3)]

    @test label_components!(zeros(UInt8, 240), trues(240); dims=()) == 1:240
    @test_throws ErrorException("labels exhausted, use a larger integer type") label_components!(
        zeros(UInt8, 260), trues(260); dims=()
    )
    @test label_components!(zeros(UInt16, 260), trues(260); dims=()) == 1:260
end
