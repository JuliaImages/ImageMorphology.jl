@testset "Connected components" begin
    @testset "label_components" begin
        @testset "interface" begin
            for T in (Bool, N0f8, Float64, Gray{N0f8})
                A = rand(T, 11, 11)
                @test label_components(A) == label_components(A, strel_diamond(A))
                @test label_components(A; r=2, dims=1) == label_components(A, strel_diamond(A, 1; r=2))

                label = similar(A, Int)
                label_components!(label, A)
                @test label_components(A) == label
                label_components!(label, A; r=2, dims=1)
                @test label_components(A; r=2, dims=1) == label
                label_components!(label, A, strel_box((3, 3)))
                @test label_components(A, strel_box((3, 3))) == label

                Ac = complement.(A)
                @test label_components(Ac; bkg=oneunit(T)) == label_components(A; bkg=zero(T))
            end

            A = rand(Bool, 11, 11)
            se = strel(CartesianIndex, centered(Bool[1 1 1; 1 1 0; 0 0 0]))
            label = @suppress_err label_components(A, se) # deprecated usage in v0.4
            @test label == label_components(A, strel_box((3, 3)))

            msg = "Non-symmetric structuring element is not supported yet"
            se = centered(Bool[1 1 1; 1 1 1; 1 1 0])
            @test_throws ArgumentError(msg) label_components(rand(Bool, 5, 5), se)

            msg = "axes of input and output must match, got (Base.OneTo(11), Base.OneTo(11)) and (-1:1, -1:1). The second argument seems to be a structuring element, it is expected to be the input array."
            @test_throws DimensionMismatch(msg) label_components!(rand(Bool, 11, 11), strel_diamond((3, 3)))

            msg = "labels exhausted, use a larger integer type"
            @test_throws ErrorException(msg) label_components!(zeros(Bool, 4, 4), rand(4, 4))
        end

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
        connectivity = centered([
            false true  false
            true  false true
            false true  false
        ])
        @test label_components(A, connectivity) == lbltarget
        connectivity = centered(trues(3, 3))
        lbltarget2 = [
            1 1 0 1
            1 0 1 1
        ]
        @test label_components(A, connectivity) == lbltarget2
    end

    @testset "traits" begin
        A = [
            true  true  false true
            true  false true  true
        ]
        label = label_components(A)
        @test component_boxes(label) ==
            Vector{Tuple}[[(1, 2), (2, 3)], [(1, 1), (2, 2)], [(1, 3), (2, 4)]]
        @test component_lengths(label) == [2, 3, 3]
        @test component_indices(label) == Array{Int}[[4, 5], [1, 2, 3], [6, 7, 8]]
        @test component_subscripts(label) ==
            Array{Tuple}[[(2, 2), (1, 3)], [(1, 1), (2, 1), (1, 2)], [(2, 3), (1, 4), (2, 4)]]
        @test @inferred(component_centroids(label)) ==
            Tuple[(1.5, 2.5), (4 / 3, 4 / 3), (5 / 3, 11 / 3)]

        @test label_components!(zeros(UInt8, 240), trues(240); dims=()) == 1:240
        @test_throws ErrorException("labels exhausted, use a larger integer type") label_components!(
            zeros(UInt8, 260), trues(260); dims=()
        )
        @test label_components!(zeros(UInt16, 260), trues(260); dims=()) == 1:260
    end

end
