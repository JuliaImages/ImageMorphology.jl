@testset "strel" begin
    @testset "offset to mask" begin
        se = [CartesianIndex(-1), CartesianIndex(1)]
        se_mask = @inferred strel(Bool, se)
        @test axes(se_mask) == (-1:1,)
        @test se_mask == centered(Bool[1, 1, 1])

        se = [CartesianIndex(-1, -1), CartesianIndex(1, 1)]
        @test strel(se) === se
        se_mask = @inferred strel(Bool, se)
        @test axes(se_mask) == (-1:1, -1:1)
        @test se_mask == centered(Bool[1 0 0; 0 1 0; 0 0 1])
        # internally, Bool is converted to SEMask
        se_mask = @inferred strel(ImageMorphology.SEMask{2}(), se)
        @test se_mask == centered(Bool[1 0 0; 0 1 0; 0 0 1])

        se = CartesianIndices((-1:1, -1:1))
        se_mask = @inferred strel(Bool, se)
        @test se_mask == centered(trues((3, 3)))

        se = CartesianIndices((1:2, -1:1))
        se_mask = @inferred strel(Bool, se)
        @test axes(se_mask) == (-2:2, -1:1)
        @test se_mask == centered(Bool[0 0 0; 0 0 0; 0 1 0; 1 1 1; 1 1 1])

        se = @inferred strel(Bool, CartesianIndex{2}[])
        @test se == centered(reshape(Bool[1], (1, 1))) # v1.7: Bool[1;;]

        # invalid cases
        se = [CartesianIndex(1), CartesianIndex(1, 1)]
        @test_throws ErrorException strel(se)
    end

    @testset "mask to offset" begin
        se = centered(Bool[1, 0, 1])
        se_offsets = @inferred strel(CartesianIndex, se)
        @test se_offsets == [CartesianIndex(-1), CartesianIndex(1)]
        @test se_offsets == strel(CartesianIndex, centered(Bool[1, 1, 1]))

        se = centered(Bool[1 0 0; 0 1 0; 0 0 1])
        se_offsets = @inferred strel(CartesianIndex, se)
        @test se_offsets == [CartesianIndex(-1, -1), CartesianIndex(1, 1)]
        # internally, CartesianIndex is converted to SEOffset
        se_offsets = @inferred strel(ImageMorphology.SEOffset{2}(), se)
        @test se_offsets == [CartesianIndex(-1, -1), CartesianIndex(1, 1)]

        se = centered(trues((3, 3)))
        se_offsets = @inferred strel(CartesianIndex, se)
        @test se_offsets == filter(x -> !iszero(x), vec(CartesianIndices((-1:1, -1:1))))

        se = @inferred strel(CartesianIndex, centered(falses(3, 3)))
        @test isempty(se)

        # test deprecation
        msg = @capture_err strel(CartesianIndex, trues(3, 3))
        @test isempty(msg) || contains(msg, "connectivity mask is expected to be a centered bool array")
        # but only for 1-based array
        se = OffsetArray(trues(3, 3), -1, -1)
        msg = "`connectivity` must be symmetric bool array"
        @test_throws ArgumentError(msg) strel(CartesianIndex, se)
    end

    @testset "center point" begin
        for se in centered.([Bool[1 0 0; 0 1 0; 0 0 1], Bool[1 0 0; 0 0 0; 0 0 1]])
            # center point is always excluded in offsets
            se_offsets = @inferred strel(CartesianIndex, se)
            @test se_offsets == [CartesianIndex(-1, -1), CartesianIndex(1, 1)]
            # but included in mask
            se_mask = @inferred strel(Bool, se_offsets)
            @test se_mask == centered(Bool[1 0 0; 0 1 0; 0 0 1])
        end
    end
end

@testset "strel_chain" begin
    se1 = centered(Bool[0 1 0; 0 1 0; 0 1 0])
    se2 = centered(Bool[0 0 0; 1 1 1; 0 0 0])
    se = @inferred strel_chain(se1, se2)
    @test se isa ImageMorphology.SEChainArray
    @test eltype(se) == Bool
    @test axes(se) == (-1:1, -1:1)
    @test se == strel_box((3, 3))
    @test se == strel_chain([se1, se2]) == strel_chain((se1, se2))
    @test se === strel_chain(se)
    @test strel_size(se) == (3, 3)

    for f in Any[dilate, erode]
        se1 = centered(rand(Bool, 3, 3))
        se2 = centered(rand(Bool, 3, 3))
        img = rand(64, 64)
        out_pipe = f(f(img, se1), se2)
        out_se = f(img, strel_chain(se1, se2))
        R = CartesianIndices(img)[2:(end - 1), 2:(end - 1)] # inner region
        @test out_pipe[R] == out_se[R]
    end
end

@testset "strel_product" begin
    se1 = centered(Bool[1, 1, 1])
    se2 = centered(Bool[1, 1, 1])
    se = strel_product(se1, se2)
    @test_broken @inferred strel_product(se1, se2) # if necessary, fix this inference issue
    @test se isa ImageMorphology.SEChainArray
    @test se == strel_box((3, 3))

    se_list = [centered(rand(Bool, sz...)) for sz in Any[(3,), (3, 3), (3, 3, 3)]]
    se = strel_product(se_list)
    @test ndims(se) == 6
end

@testset "strel_diamond" begin
    @testset "N=1" begin
        img = rand(5)
        se = @inferred strel_diamond(img)
        @test se isa ImageMorphology.SEDiamondArray
        @test eltype(se) == Bool
        @test axes(se) == (-1:1,)
        @test se == centered(Bool[1, 1, 1])

        se = @inferred strel_diamond((5,); r=1)
        @test se == centered(Bool[0, 1, 1, 1, 0])

        se = @inferred strel_diamond((5,))
        @test se == centered(Bool[1, 1, 1, 1, 1])
    end

    @testset "N=2" begin
        img = rand(5, 5)
        se = @inferred strel_diamond(img)
        @test se isa ImageMorphology.SEDiamondArray
        @test eltype(se) == Bool
        @test axes(se) == (-1:1, -1:1)
        @test se == centered(Bool[0 1 0; 1 1 1; 0 1 0])
        @test se == strel_diamond((3, 3), (1, 2); r=1)
        @test se == strel_diamond(img, (1, 2); r=1)

        se = @inferred strel_diamond(img, (1,))
        @test se == @inferred strel_diamond(img, 1)
        @test se == centered(reshape(Bool[1, 1, 1], (3, 1)))

        se = @inferred strel_diamond((5, 5), (2, 1))
        @test se == strel_diamond((5, 5), 1:2) # not inferable

        se = @inferred strel_diamond((5, 5))
        @test se == centered(Bool[0 0 1 0 0; 0 1 1 1 0; 1 1 1 1 1; 0 1 1 1 0; 0 0 1 0 0])

        se = @inferred strel_diamond((3, 5))
        @test se == centered(Bool[0 1 1 1 0; 1 1 1 1 1; 0 1 1 1 0])

        se = @inferred strel_diamond((3, 5), (1,))
        @test se == @inferred strel_diamond((3, 5), 1)
        @test se == centered(Bool[0 0 1 0 0; 0 0 1 0 0; 0 0 1 0 0])

        se = @inferred strel_diamond((3, 5), (2,))
        @test se == centered(Bool[0 0 0 0 0; 1 1 1 1 1; 0 0 0 0 0])

        se = @inferred strel_diamond((3, 5); r=2)
        @test se == centered(Bool[0 1 1 1 0; 1 1 1 1 1; 0 1 1 1 0])
    end

    @testset "N=3" begin
        img = rand(5, 5, 5)
        se = @inferred strel_diamond(img)
        @test se isa ImageMorphology.SEDiamondArray
        @test eltype(se) == Bool
        @test axes(se) == (-1:1, -1:1, -1:1)
        @test se[:, :, -1] == se[:, :, 1] == centered(Bool[0 0 0; 0 1 0; 0 0 0])
        @test se[:, :, 0] == strel_diamond((3, 3))

        se = @inferred strel_diamond((3, 3, 3), (1, 2))
        @test se[:, :, -1] == se[:, :, 1] == centered(falses((3, 3)))
        @test se[:, :, 0] == strel_diamond((3, 3))
    end

    @testset "strel conversion" begin
        se = strel_diamond((3, 3))
        se_offsets = @inferred strel(CartesianIndex, se)
        @test se_offsets == [CartesianIndex(0, -1), CartesianIndex(-1, 0), CartesianIndex(1, 0), CartesianIndex(0, 1)]
        se_mask = @inferred strel(Bool, se)
        # not BitMatrix as SEDiamondArray provides more information of what the SE is
        @test se_mask isa ImageMorphology.SEDiamondArray
        @test se_mask === se
    end

    # edge cases
    img = rand(5, 5)
    err = ArgumentError("`axes` length should be at least 2")
    @test_throws err strel_diamond((3,), (1, 2))
    err = ArgumentError("size should be odd integers")
    @test_throws err strel_diamond((2, 3))
    err = ArgumentError("dims should be unique")
    @test_throws err strel_diamond((3, 3), (1, 1))
    err = ArgumentError("all `dims` values should be less than or equal to 2")
    @test_throws err strel_diamond((3, 3), (5,))
end

@testset "strel_box" begin
    @testset "N=1" begin
        img = rand(5)
        se = @inferred strel_box(img)
        @test se isa ImageMorphology.SEBoxArray
        @test eltype(se) == Bool
        @test se == centered(Bool[1, 1, 1])

        se = @inferred strel_box((5,))
        @test se == centered(Bool[1, 1, 1, 1, 1])

        se = @inferred strel_box((5,); r=1)
        @test se == centered(Bool[0, 1, 1, 1, 0])
    end

    @testset "N=2" begin
        img = rand(5, 5)
        se = @inferred strel_box(img)
        @test se isa ImageMorphology.SEBoxArray
        @test eltype(se) == Bool
        @test se == centered(Bool[1 1 1; 1 1 1; 1 1 1])
        @test se == strel_box(img; r=1)

        se = @inferred strel_box(img, (1,))
        @test se == @inferred strel_box((3, 1), (1,))
        @test se == @inferred strel_box((3, 1), 1)

        se = @inferred strel_box(img; r=2)
        @test se == centered(trues((5, 5)))

        se = @inferred strel_box(img; r=(1, 2))
        @test se == centered(trues((3, 5)))

        se = @inferred strel_box((5, 5), (2, 1))
        @test se == @inferred strel_box((5, 5), 1:2)

        se = @inferred strel_box((3, 5))
        @test se == centered(Bool[1 1 1 1 1; 1 1 1 1 1; 1 1 1 1 1])

        se = @inferred strel_box((3, 5); r=1)
        @test se == centered(Bool[0 1 1 1 0; 0 1 1 1 0; 0 1 1 1 0])

        se = @inferred strel_box((3, 5); r=(1, 0))
        @test se == @inferred strel_box((3, 5), (1,))
        @test se == @inferred strel_box((3, 5), 1)
        @test se == centered(Bool[0 0 1 0 0; 0 0 1 0 0; 0 0 1 0 0])

        se = @inferred strel_box((3, 5); r=(0, 1))
        @test se == centered(Bool[0 0 0 0 0; 0 1 1 1 0; 0 0 0 0 0])

        se = @inferred strel_box((3, 5), (2,))
        @test se == centered(Bool[0 0 0 0 0; 1 1 1 1 1; 0 0 0 0 0])

        se = @inferred strel_box((3, 5); r=2)
        @test se == centered(Bool[1 1 1 1 1; 1 1 1 1 1; 1 1 1 1 1])
    end

    @testset "N=3" begin
        img = rand(5, 5, 5)
        se = @inferred strel_box(img)
        @test se isa ImageMorphology.SEBoxArray
        @test eltype(se) == Bool
        @test se[:, :, -1] == se[:, :, 1] == centered(trues((3, 3)))
        @test se[:, :, 0] == strel_box((3, 3))

        se = @inferred strel_box((3, 3, 3), (1, 2))
        @test axes(se) == (-1:1, -1:1, -1:1)
        @test se[:, :, 0] == strel_box((3, 3))
    end

    @testset "strel conversion" begin
        se = strel_box((3, 3))
        se_offsets = @inferred strel(CartesianIndex, se)
        @test se_offsets == filter(i -> !iszero(i), vec(CartesianIndices((-1:1, -1:1))))
        se_mask = @inferred strel(Bool, se)
        # not BitMatrix as SEBoxArray provides more information of what the SE is
        @test se_mask isa ImageMorphology.SEBoxArray
        @test se_mask === se
    end
end

@testset "strel_split" begin
    se = strel_diamond((3, 3))
    upper, lower = ImageMorphology.strel_split(se)
    @test upper == centered(Bool[0 1 0; 1 1 0; 0 0 0])
    @test lower == centered(Bool[0 0 0; 0 1 1; 0 1 0])
    se = strel(CartesianIndex, se)
    upper, lower = ImageMorphology.strel_split(se)
    @test upper == [CartesianIndex(0, -1), CartesianIndex(-1, 0)]
    @test lower == [CartesianIndex(1, 0), CartesianIndex(0, 1)]

    # non-symmetric is not allowed
    se = collect(strel_diamond((3, 3)))
    se[1, 1] = true
    msg = "structuring element must be symmetric with respect to its center"
    @test_throws ArgumentError(msg) ImageMorphology.strel_split(se)
end

@testset "strel_type" begin
    @test strel_type(strel_diamond((3, 3))) isa ImageMorphology.SEDiamond
    @test strel_type(strel_box((3, 3))) isa ImageMorphology.SEBox
    @test strel_type([CartesianIndex(1, 2)]) isa ImageMorphology.SEOffset
    @test strel_type(CartesianIndices((-1:1, -1:1))) isa ImageMorphology.SEOffset
    @test strel_type(trues(3, 3)) isa ImageMorphology.SEMask

    msg = "invalid structuring element data type: $(Vector{CartesianIndex})"
    @test_throws ErrorException(msg) strel_type(CartesianIndex[])
end

@testset "strel_size" begin
    se = strel_diamond((5, 5), (1,))
    @test strel_size(se) == strel_size(centered(collect(se))) == (5, 1)

    se = strel_box((5, 5), (2,))
    @test strel_size(se) == strel_size(centered(collect(se))) == (1, 5)

    se = [CartesianIndex(-2, -2), CartesianIndex(1, 1)]
    @test strel_size(se) == (5, 5)
end

@testset "centered" begin
    se = strel_diamond((3, 3))
    @test centered(se) === se
    se = strel_box((3, 3))
    @test centered(se) === se
end
