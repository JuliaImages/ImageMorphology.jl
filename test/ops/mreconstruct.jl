ref_mreconstruct(::typeof(dilate), args...; kwargs...) = ref_mreconstruct((dilate, min), args...; kwargs...)
ref_mreconstruct(::typeof(erode), args...; kwargs...) = ref_mreconstruct((erode, max), args...; kwargs...)
function ref_mreconstruct((op, select), marker, mask; dims=coords_spatial(mask))
    return ref_mreconstruct((op, select), marker, mask, strel_box(marker, dims))
end
function ref_mreconstruct((op, select), marker, mask, se)
    # maximum iteration limit to spread top-left pixel to entire image
    max_iter = maximum(size(mask))
    out = select.(mask, marker)
    has_changed = true
    n = 1
    while has_changed && n < max_iter
        old = out
        out = select.(op(out, se), mask)
        has_changed = !(old == out) # must be == instead of ≈
        n += 1
    end
    return out
end

@testset "mreconstruct" begin
    @testset "API" begin
        mask, marker = rand(11, 11), rand(11, 11)
        out1 = mreconstruct(dilate, marker, mask)
        out2 = similar(mask)
        out2 = mreconstruct!(dilate, out2, marker, mask)
        @test out1 == out2
        @test axes(out1) == (1:11, 1:11)

        # ensure the eltype can be mixed
        mask, marker = rand(Gray{Float64}, 11, 11), rand(N0f8, 11, 11)
        out1 = mreconstruct(dilate, marker, mask)
        @test eltype(out1) == Gray{Float64}
        out2 = mreconstruct(dilate, mask, marker)
        @test out1 != out2 # the result changes if swapping the order
        @test eltype(out2) == Gray{Float64}

        se = strel_diamond((7, 7))
        msg = @capture_err mreconstruct(dilate, marker, mask, se)
        @test occursin("structuring element with half-size larger than 1 is invalid, only the center 3×3 values are used", msg)
        ref = mreconstruct(dilate, marker, mask)
        @test ref == @suppress_err mreconstruct(dilate, marker, mask, se)

        msg = "operation `max` is not supported for `mreconstruct`"
        @test_throws ErrorException(msg) mreconstruct(max, marker, mask)
        msg = "`marker` and `mask` should have the same axes"
        @test_throws DimensionMismatch(msg) mreconstruct(dilate, rand(10), rand(3))
        msg = "the input structuring element is not for 1 dimensional array, instead it is for 2 dimensional array"
        @test_throws DimensionMismatch(msg) mreconstruct(dilate, rand(10), rand(10), strel_box((3, 3)))
        @test_throws MethodError mreconstruct(dilate, rand(10), rand(10), strel_box((3,)); dims=1)
    end

    @testset "numeric" begin
        # 1D
        marker = [0, 0, 0, 3, 0]
        mask = [0, 0, 6, 5, 6]
        ref_out = [0, 0, 3, 3, 3]
        @test ref_out == ref_mreconstruct(dilate, marker, mask)
        out = mreconstruct(dilate, marker, mask)
        @test ref_out == out

        # 2D -- marker > mask for non-zero values
        marker = [0 0 0 0 0; 0 9 0 0 0; 0 0 0 0 0; 0 0 0 5 0; 0 0 0 0 0; 0 9 0 0 0]
        mask = [9 0 0 0 0; 0 8 7 1 0; 0 9 0 4 0; 0 0 0 4 0; 0 0 6 5 6; 0 0 9 8 9]
        ref_out = [8 0 0 0 0; 0 8 7 1 0; 0 8 0 4 0; 0 0 0 4 0; 0 0 4 4 4; 0 0 4 4 4]
        out = mreconstruct(dilate, marker, mask)
        @test out == ref_out
        @test ref_out == ref_mreconstruct(dilate, marker, mask)
        ref_out = [0 0 0 0 0; 0 8 7 1 0; 0 8 0 4 0; 0 0 0 4 0; 0 0 4 4 4; 0 0 4 4 4]
        out = mreconstruct(dilate, marker, mask, strel_diamond(marker))
        @test out == ref_out
        @test ref_out == ref_mreconstruct(dilate, marker, mask, strel_diamond(marker))

        # 2D -- marker > mask for non-zero values doesn't hold
        marker = [0 0 0 0 0; 0 9 0 0 0; 0 0 0 0 0; 0 0 0 3 0; 0 0 0 0 0; 0 9 0 0 0]
        mask = [0 0 0 0 0; 0 8 7 1 0; 0 9 0 4 0; 0 0 0 4 0; 0 0 6 5 6; 0 0 9 8 9]
        ref_out = [0 0 0 0 0; 0 8 7 1 0; 0 8 0 4 0; 0 0 0 4 0; 0 0 4 4 4; 0 0 4 4 4]
        out = mreconstruct(dilate, marker, mask)
        @test out == ref_out

        # 3D
        marker = zeros(5, 6, 3)
        marker[3, 4, 2] = 3

        mask = zeros(5, 6, 3)
        mask[3, 4, 1] = 4
        mask[2, 3:5, 2] = [5, 5, 5]
        mask[3, 3:5, 2] = [4, 4, 4]
        mask[4, 3:5, 2] = [4, 4, 4]
        mask[3, 4, 2] = 6
        mask[3, 4, 3] = 4

        ref_out = zeros(5, 6, 3)
        ref_out[3, 4, 1] = 3
        ref_out[2, 3:5, 2] = [3, 3, 3]
        ref_out[3, 3:5, 2] = [3, 3, 3]
        ref_out[4, 3:5, 2] = [3, 3, 3]
        ref_out[3, 4, 3] = 3
        out = mreconstruct(dilate, marker, mask)
        @test out == ref_out
        @test ref_out == ref_mreconstruct(dilate, marker, mask)

        # extensive numeric test against different types
        for T in Any[Bool, Int, N0f8, Gray{N0f8}, Gray{Float64}, Float64]
            for N in (1, 2, 3)
                sz = ntuple(_ -> 64, N)
                mask = T == Int ? rand(1:10, sz...) : rand(T, sz...)
                dims_list = ntuple(i -> ntuple(identity, i), N)
                for dims in dims_list
                    marker = T == Int ? rand(1:10, sz...) : rand(T, sz...)
                    for op in Any[erode, dilate]
                        ref_out = ref_mreconstruct(op, marker, mask; dims)
                        out = mreconstruct(op, marker, mask; dims)
                        @test out == ref_out

                        outc = mreconstruct(op, centered(marker), centered(mask); dims)
                        @test outc == centered(out)
                    end

                    # duality
                    mask = T == Int ? rand(1:10, sz...) : rand(T, sz...)
                    marker = T == Int ? rand(1:10, sz...) : rand(T, sz...)
                    out1 = mreconstruct(dilate, marker, mask)
                    out2 = mreconstruct(erode, complement.(marker), complement.(mask))
                    @test out1 == complement.(out2)
                end
            end
        end
    end
end
