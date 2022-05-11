@testset "extreme_filter" begin
    @testset "numerical" begin
        # Bool
        img = fill(false, 5, 5); img[3, 3] = 1
        ref = collect(strel_diamond((5, 5); r=1))

        img_d = @inferred extreme_filter(max, img)
        @test img_d == ref
        @test ref == .!extreme_filter(min, .!img) # dual property

        # Int
        A = Int[4 6 5 3 4; 8 6 9 4 8; 7 8 4 9 6; 6 2 2 1 7; 1 6 5 2 6]
        ref = [8 6 9 5 8; 8 9 9 9 8; 8 8 9 9 9; 7 8 5 9 7; 6 6 6 6 7]
        Ad = @inferred extreme_filter(max, A)
        @test Ad == ref

        # Gray{Float32}
        img = Gray{Float32}.(A./9)
        img_d = @inferred extreme_filter(max, img)
        @test ref == gray.(img_d).*9
    end

    @testset "strel" begin
        # ensure it supports both strel representations
        A = rand(32, 32)
        se_mask = centered(collect(strel_diamond(A)))
        out = extreme_filter(max, A, se_mask)
        se_offsets = strel(CartesianIndex, se)
        @test out == extreme_filter(max, A, se_offsets)
    end

    @testset "diamond" begin
        # ensure the optimized implementation work equivalently to the generic fallback implementation
        for N in (1, 2, 3)
            sz = ntuple(_->32, N)
            img = rand(sz...)
            for r in (1, 3)
                dims_list = ntuple(i->ntuple(identity, i), N)
                for dims in dims_list
                    se = centered(collect(strel_diamond(img, dims)))
                    ref = extreme_filter(max, img, se) # this hits the generic path

                    out = extreme_filter(max, img; dims=dims) # this hits the optimized path
                    @test out == ref
                end
            end
        end
    end
end
