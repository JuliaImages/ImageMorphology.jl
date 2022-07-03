@testset "overbuild" begin
    # heavy numerical tests are put in the `mreconstruct` testset
    marker, mask = rand(32, 32), rand(32, 32)
    @test overbuild(marker, mask) ==
        mreconstruct(erode, marker, mask)
    @test overbuild(marker, mask; dims=1) ==
        mreconstruct(erode, marker, mask; dims=1)
    @test overbuild(marker, mask, strel_diamond(marker)) ==
        mreconstruct(erode, marker, mask, strel_diamond(marker))
    @test overbuild!(similar(marker), marker, mask) ==
        mreconstruct!(erode, similar(marker), marker, mask)
    @test overbuild!(similar(marker), marker, mask; dims=1) ==
        mreconstruct!(erode, similar(marker), marker, mask; dims=1)
    @test overbuild!(similar(marker), marker, mask, strel_diamond(marker)) ==
        mreconstruct!(erode, similar(marker), marker, mask, strel_diamond(marker))
end
