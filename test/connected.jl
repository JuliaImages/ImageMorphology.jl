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

    @testset "component_boxes" begin
        A = rand(0:5, 11, 11)
        label = label_components(A, strel_box((5, 5)))
        boxes = component_boxes(label)
        @test eltype(boxes) <: CartesianIndices
        @test axes(boxes) == (0:maximum(label),)

        # the bounding box is the smallest region that contains the label value
        @test all(eachindex(boxes)) do i
            n1 = count(isequal(i), label[boxes[i]])
            n2 = count(isequal(i), label)
            n1 == n2
        end

        A = [2 2 2 2 2; 1 1 1 0 1; 1 0 2 1 1; 1 1 2 2 2; 1 0 2 2 2]
        label = label_components(A)
        boxes = component_boxes(label)
        @test boxes == OffsetArray(
            [
                CartesianIndices((2:5, 2:4)),
                CartesianIndices((1:1, 1:5)),
                CartesianIndices((2:5, 1:3)),
                CartesianIndices((3:5, 3:5)),
                CartesianIndices((2:3, 4:5)),
            ], -1)
        @test A[boxes[1]] == [2 2 2 2 2]
        @test A[boxes[4]] == [0 1; 1 1]

        label = zeros(Int, 5, 5)
        label[1] = -1
        msg = "The input labeled array should contain background label `0` as the minimum value"
        @test_throws ArgumentError(msg) component_boxes(label)

        boxes = component_boxes(ones(Int, 5, 5))
        @test length(boxes) == 1
    end

    @testset "component_indices" begin
        A = [2 2 2 2 2; 1 1 1 0 1; 1 0 2 1 1; 1 1 2 2 2; 1 0 2 2 2]
        label = label_components(A)
        indices = component_indices(label)
        @test IndexStyle(label) == IndexLinear()
        @test eltype(indices) == Vector{Int}
        @test axes(indices) == (0:maximum(label),)
        @test indices == OffsetVector([
                [8, 10, 17], [1, 6, 11, 16, 21], [2, 3, 4, 5, 7, 9, 12],
                [13, 14, 15, 19, 20, 24, 25], [18, 22, 23],
            ], -1)

        @test indices == component_indices(Int, label) == component_indices(IndexLinear(), label)

        indices = component_indices(CartesianIndex, label)
        @test eltype(indices) == Vector{CartesianIndex{ndims(A)}}
        @test axes(indices) == (0:maximum(label),)
        @test indices == [CartesianIndices(label)[R] for R in component_indices(Int, label)]
        @test indices == component_indices(IndexCartesian(), label)

        plabel = PaddedView(0, label, (0:6, 0:6))
        indices = component_indices(CartesianIndex, label)
        pindices = component_indices(plabel)
        @test IndexStyle(plabel) == IndexCartesian()
        @test eltype(pindices) == Vector{CartesianIndex{ndims(A)}}
        @test pindices[1:end] == indices[1:end] # they only differ in background indices

        label = zeros(Int, 5, 5)
        label[1] = -1
        msg = "The input labeled array should contain background label `0` as the minimum value"
        @test_throws ArgumentError(msg) component_indices(label)

        indices = component_indices(ones(Int, 5, 5))
        @test axes(indices, 1) == 0:1
        @test isempty(indices[0])
    end

    @testset "component_lengths" begin
        A = [2 2 2 2 2; 1 1 1 0 1; 1 0 2 1 1; 1 1 2 2 2; 1 0 2 2 2]
        label = label_components(A)
        counts = component_lengths(label)
        indices = component_indices(label)
        @test counts == length.(indices)
        @test axes(counts) == (0:maximum(label),)

        label = zeros(Int, 5, 5)
        label[1] = -1
        msg = "The input labeled array should contain background label `0` as the minimum value"
        @test_throws ArgumentError(msg) component_lengths(label)

        counts = component_lengths(ones(Int, 5, 5))
        @test axes(counts, 1) == 0:1
        @test counts[0] == 0
    end

    @testset "component_centroids" begin
        A = [2 2 2 2 2; 1 1 1 0 1; 1 0 2 1 1; 1 1 2 2 2; 1 0 2 2 2]
        label = label_components(A)
        centers = component_centroids(label)
        indices = component_indices(CartesianIndex, label)
        @test centers == [sum(X).I ./ length(X) for X in indices]
        @test axes(centers) == (0:maximum(label),)

        label = zeros(Int, 5, 5)
        label[1] = -1
        msg = "The input labeled array should contain background label `0` as the minimum value"
        @test_throws ArgumentError(msg) component_centroids(label)

        centers = component_centroids(ones(Int, 5, 5))
        @test axes(centers, 1) == 1:1
    end

    @testset "deprecations" begin
        A = [
            true  true  false true
            true  false true  true
        ]
        label = label_components(A)
        indices = @suppress_err component_subscripts(label)
        @test indices ==
            Array{Tuple}[[(2, 2), (1, 3)], [(1, 1), (2, 1), (1, 2)], [(2, 3), (1, 4), (2, 4)]]
    end

end
