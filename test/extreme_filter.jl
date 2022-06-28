@testset "extreme_filter" begin
    @testset "interface" begin
        img = rand(1:10, 7, 7)

        out = extreme_filter(max, img)
        @test out == extreme_filter(max, img; dims=(1, 2), r=1)

        se = strel_box(img, (1, 2); r=(1, 1))
        @test out == extreme_filter(max, img, se)

        se = strel_box(img, (1, 2); r=(1, 3))
        @test extreme_filter(max, img; dims=(1, 2), r=(1, 3)) == extreme_filter(max, img, se)

        out = similar(img)
        rst = extreme_filter!(max, out, img)
        @test rst === out
        @test out == extreme_filter!(max, out, img)

        out = similar(img)
        rst = extreme_filter!(max, out, img; dims=(1, 2), r=1)
        @test rst === out
        se = strel_box(img, (1, 2); r=(1, 1))
        ref_out = similar(img)
        extreme_filter!(max, ref_out, img, se)
        @test ref_out == out

        out = similar(img)
        rst = extreme_filter!(max, out, img; dims=(1, 2), r=(1, 3))
        @test rst === out
        se = strel_box(img, (1, 2); r=(1, 3))
        ref_out = similar(img)
        extreme_filter!(max, ref_out, img, se)
        @test ref_out == out
    end

    @testset "numerical" begin
        # Bool
        img = fill(false, 5, 5)
        img[3, 3] = 1
        ref = collect(strel_box((5, 5); r=1))

        img_d = @inferred extreme_filter(max, img)
        @test img_d == ref
        @test ref == .!extreme_filter(min, .!img) # dual property

        # Int
        A = Int[4 6 5 3 4; 8 6 9 4 8; 7 8 4 9 6; 6 2 2 1 7; 1 6 5 2 6]
        ref = [8 6 9 5 8; 8 9 9 9 8; 8 8 9 9 9; 7 8 5 9 7; 6 6 6 6 7]
        Ad = @inferred extreme_filter(max, A, strel_diamond(A))
        @test Ad == ref
        ref = [8 9 9 9 8; 8 9 9 9 9; 8 9 9 9 9; 8 8 9 9 9; 6 6 6 7 7]
        Ad = @inferred extreme_filter(max, A)
        @test Ad == ref

        # Gray{Float32}
        img = Gray{Float32}.(A ./ 9)
        img_d = @inferred extreme_filter(max, img)
        @test ref == gray.(img_d) .* 9
    end

    @testset "strel" begin
        # ensure it supports both strel representations
        A = rand(32, 32)
        se_mask = centered(collect(strel_diamond(A)))
        out = extreme_filter(max, A, se_mask)
        se_offsets = strel(CartesianIndex, se_mask)
        @test out == extreme_filter(max, A, se_offsets)
    end

    @testset "offset arrays" begin
        img = centered(rand(1:10, 32, 32))
        out = extreme_filter(max, img)
        @test out == centered(extreme_filter(max, OffsetArrays.no_offset_view(img)))
    end

    @testset "multi channel image" begin
        # supporting color image requires some reduced ordering function
        img = rand(RGB, 32, 32)
        msg = "function `max` is not a well-defined select function on type `RGB{Float64}`: does `f(x::T, y::T)` work as expected?"
        @test_throws ArgumentError(msg) extreme_filter(max, img)

        # the test functions below don't make much sense in practice, but it's good to test them
        _to_number(c::AbstractRGB) = red(c)
        _to_number(x::Number) = x
        _select(x, y) = _to_number(x) > _to_number(y) ? x : y
        out = extreme_filter(_select, img)
        @test issubset(Set(unique(out)), Set(unique(img)))

        _select_new(x, y) = max(_to_number(x), _to_number(y))
        out = extreme_filter(_select_new, img)
        @test eltype(out) == eltype(img)
        @test sum(abs, channelview(RGB.(Gray.(out))) - channelview(out)) < 1e-4
    end

    # Ensure our various implementations of extreme_filter are equivalent by testing against
    # the generic implementation
    for (se_gen_func, se_name) in Any[(strel_diamond, "diamond"), (strel_box, "box")]
        test_name = "optimization: $se_name"
        @testset "$test_name" begin
            # ensure the optimized implementation work equivalently to the generic fallback implementation
            for T in Any[Bool, Int, N0f8, Gray{N0f8}, Gray{Float64}, Float64]
                for N in (1, 2, 3)
                    sz = ntuple(_ -> 32, N)
                    img = T == Int ? rand(1:10, sz...) : rand(T, sz...)
                    for r in (1, 3)
                        dims_list = ntuple(i -> ntuple(identity, i), N)
                        for dims in dims_list
                            se = se_gen_func(ntuple(_ -> 2r + 1, N), dims)
                            ref = ImageMorphology._extreme_filter_generic!(max, similar(img), img, se)
                            out = extreme_filter(max, img, se)
                            @test out == ref

                            # `maybe_floattype` is used to pre-allocate the output for
                            # higher level operators (e.g., tophat) to avoid overflow
                            # behavior. We need to ensure it works correctly.
                            # https://github.com/JuliaImages/ImageMorphology.jl/issues/104
                            out = similar(img, ImageMorphology.maybe_floattype(T))
                            extreme_filter!(max, out, img, se)
                            @test out == ref

                            imgc = centered(img)
                            outc = extreme_filter(max, imgc, se)
                            @test outc == centered(out)
                        end
                    end
                end
            end
        end
    end

    @testset "optimization: bool" begin
        # ensure the optimized implementation work equivalently to the generic fallback implementation
        for N in (1, 2, 3)
            sz = ntuple(_ -> 32, N)
            img = rand(Bool, sz...)
            for r in (1, 3)
                dims_list = ntuple(i -> ntuple(identity, i), N)
                for dims in dims_list, select in (max, min)
                    se = strel_diamond(ntuple(_ -> 2r + 1, N), dims)
                    ref = ImageMorphology._extreme_filter_generic!(select, similar(img), img, se)
                    out = extreme_filter(select, img, se) # SEDiamond method
                    @test out == ref

                    se = centered(collect(se))
                    out = extreme_filter(select, img, se) # MorphologySE method
                    @test out == ref
                end
            end
        end
    end
end
