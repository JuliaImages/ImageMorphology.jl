@testset "underbuild" begin
    # heavy numerical tests are put in the `mreconstruct` testset
    marker, mask = rand(32, 32), rand(32, 32)
    @test underbuild(marker, mask) ==
        mreconstruct(dilate, marker, mask)
    @test underbuild(marker, mask; dims=1) ==
        mreconstruct(dilate, marker, mask; dims=1)
    @test underbuild(marker, mask, strel_diamond(marker)) ==
        mreconstruct(dilate, marker, mask, strel_diamond(marker))
    @test underbuild!(similar(marker), marker, mask) ==
        mreconstruct!(dilate, similar(marker), marker, mask)
    @test underbuild!(similar(marker), marker, mask; dims=1) ==
        mreconstruct!(dilate, similar(marker), marker, mask; dims=1)
    @test underbuild!(similar(marker), marker, mask, strel_diamond(marker)) ==
        mreconstruct!(dilate, similar(marker), marker, mask, strel_diamond(marker))
end
