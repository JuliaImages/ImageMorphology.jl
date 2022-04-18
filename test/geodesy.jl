@testset "underbuild" begin
    #1D
    marker = [0 0 0 3 0]
    mask = [0 0 6 5 6]
    expected_connectivity_1D = [0 0 3 3 3]
    connectivity_1D = trues(1, 3)
    output_1D = underbuild(marker, mask, connectivity_1D)
    @test output_1D == expected_connectivity_1D

    #2D
    marker = ([0 0 0 0 0
        0 9 0 0 0
        0 0 0 0 0
        0 0 0 3 0
        0 0 0 0 0
        0 9 0 0 0])

    mask = ([0 0 0 0 0
        0 8 7 1 0
        0 9 0 4 0
        0 0 0 4 0
        0 0 6 5 6
        0 0 9 8 9])

    expected_connectivity_2D8 = ([0 0 0 0 0
        0 8 7 1 0
        0 8 0 4 0
        0 0 0 4 0
        0 0 4 4 4
        0 0 4 4 4])
    connectivity_2D8 = trues(3, 3)
    output_2D8 = underbuild(marker, mask, connectivity_2D8)
    @test output_2D8 == expected_connectivity_2D8

    expected_connectivity_2D4 = ([0 0 0 0 0
        0 8 7 1 0
        0 8 0 3 0
        0 0 0 3 0
        0 0 3 3 3
        0 0 3 3 3])
    connectivity_2D4 = [false true false
        true false true
        false true false]
    output_2D4 = underbuild(marker, mask, connectivity_2D4)
    @test output_2D4 == expected_connectivity_2D4

    #3D
    marker = zeros(5, 6, 3)
    marker[3, 4, 2] = 3

    mask = zeros(5, 6, 3)
    mask[3, 4, 1] = 4
    mask[2, 3:5, 2] = [5, 5, 5]
    mask[3, 3:5, 2] = [4, 4, 4]
    mask[4, 3:5, 2] = [4, 4, 4]
    mask[3, 4, 2] = 6
    mask[3, 4, 3] = 4

    expected_connectivity_3D26 = zeros(5, 6, 3)
    expected_connectivity_3D26[3, 4, 1] = 3
    expected_connectivity_3D26[2, 3:5, 2] = [3, 3, 3]
    expected_connectivity_3D26[3, 3:5, 2] = [3, 3, 3]
    expected_connectivity_3D26[4, 3:5, 2] = [3, 3, 3]
    expected_connectivity_3D26[3, 4, 3] = 3
    connectivity_3D26 = trues(3, 3, 3)
    output_3D26 = underbuild(marker, mask, connectivity_3D26)
    @test output_3D26 == expected_connectivity_3D26

    #specilization binary cases
    marker = ([
        0 0 0 0 0 0
        0 0 0 0 0 0
        0 0 0 0 0 0
        0 0 0 0 0 0
        0 0 0 0 0 0
        0 0 0 0 0 1])
    marker = convert(Array{Bool}, marker)
    mask = ([
        0 0 0 0 0 0
        0 1 1 0 0 0
        0 1 1 0 0 0
        0 1 1 0 0 0
        0 0 0 1 1 1
        0 0 0 1 1 1])
    mask = convert(Array{Bool}, mask)

    expected_connectivity_2D8 = ([
        0 0 0 0 0 0
        0 1 1 0 0 0
        0 1 1 0 0 0
        0 1 1 0 0 0
        0 0 0 1 1 1
        0 0 0 1 1 1])
    expected_connectivity_2D8 = convert(Array{Bool}, expected_connectivity_2D8)
    connectivity_2D8 = trues(3, 3)
    output_2D8 = underbuild(marker, mask, connectivity_2D8)
    @test output_2D8 == expected_connectivity_2D8

    expected_connectivity_2D4 = ([
        0 0 0 0 0 0
        0 0 0 0 0 0
        0 0 0 0 0 0
        0 0 0 0 0 0
        0 0 0 1 1 1
        0 0 0 1 1 1])
    expected_connectivity_2D4 = convert(Array{Bool}, expected_connectivity_2D4)

    connectivity_2D4 = [false true false
        true false true
        false true false]
    output_2D4 = underbuild(marker, mask, connectivity_2D4)
    @test output_2D4 == expected_connectivity_2D4
end

@testset "overbuild" begin
    marker = ([9 9 9 9 9
        9 8 9 9 9
        9 9 9 9 9
        9 9 9 3 9
        9 9 9 9 9
        9 8 9 9 9])

    mask = ([9 9 9 9 9
        9 2 3 8 9
        9 1 9 6 9
        9 9 9 6 9
        9 9 4 5 4
        9 9 1 2 1])

    expected_connectivity_2D8 = ([9 9 9 9 9
        9 6 6 8 9
        9 6 9 6 9
        9 9 9 6 9
        9 9 6 6 6
        9 9 6 6 6])
    connectivity_2D8 = trues(3, 3)
    output_2D8 = overbuild(marker, mask, connectivity_2D8)
    @test output_2D8 == expected_connectivity_2D8

    expected_connectivity_2D4 = ([9 9 9 9 9
        9 8 8 8 9
        9 8 9 6 9
        9 9 9 6 9
        9 9 6 6 6
        9 9 6 6 6])
    connectivity_2D4 = [false true false
        true false true
        false true false]
    output_2D4 = overbuild(marker, mask, connectivity_2D4)
    @test output_2D4 == expected_connectivity_2D4

    #specilization binary cases
    marker = ([
        1 1 1 1 1 1
        1 1 1 1 1 1
        1 1 1 1 1 1
        1 1 1 1 1 1
        1 1 1 1 1 1
        1 1 1 1 1 0])
    marker = convert(Array{Bool}, marker)
    mask = ([
        1 1 1 1 1 1
        1 0 0 1 1 1
        1 0 0 1 1 1
        1 0 0 1 1 1
        1 1 1 0 0 0
        1 1 1 0 0 0])
    mask = convert(Array{Bool}, mask)

    expected_connectivity_2D8 = ([
        1 1 1 1 1 1
        1 0 0 1 1 1
        1 0 0 1 1 1
        1 0 0 1 1 1
        1 1 1 0 0 0
        1 1 1 0 0 0])
    expected_connectivity_2D8 = convert(Array{Bool}, expected_connectivity_2D8)
    connectivity_2D8 = trues(3, 3)
    output_2D8 = overbuild(marker, mask, connectivity_2D8)
    @test output_2D8 == expected_connectivity_2D8

    expected_connectivity_2D4 = ([
        1 1 1 1 1 1
        1 1 1 1 1 1
        1 1 1 1 1 1
        1 1 1 1 1 1
        1 1 1 0 0 0
        1 1 1 0 0 0])
    expected_connectivity_2D4 = convert(Array{Bool}, expected_connectivity_2D4)

    connectivity_2D4 = [false true false
        true false true
        false true false]
    output_2D4 = overbuild(marker, mask, connectivity_2D4)
    @test output_2D4 == expected_connectivity_2D4
end

@testset "hmaxima" begin
    image = ([
        3 3 3 3 3 3 3 3 3 3
        3 4 4 4 3 3 4 4 4 3
        3 4 5 4 3 3 4 9 4 3
        3 4 4 4 3 3 4 4 4 3
        3 3 3 3 3 3 3 3 3 3])

    expected_connectivity_2D8 = ([
        3 3 3 3 3 3 3 3 3 3
        3 3 3 3 3 3 4 4 4 3
        3 3 3 3 3 3 4 6 4 3
        3 3 3 3 3 3 4 4 4 3
        3 3 3 3 3 3 3 3 3 3])
    connectivity_2D8 = trues(3, 3)
    output_2D8 = hmaxima(image, connectivity_2D8, (3))
    @test output_2D8 == expected_connectivity_2D8
end

@testset "hminima" begin
    image = ([
        8 8 8 8 8 8 8 8 8 8
        8 7 7 7 8 8 7 7 7 8
        8 7 6 7 8 8 7 3 7 8
        8 7 7 7 8 8 7 7 7 8
        8 8 8 8 8 8 8 8 8 8])

    expected_connectivity_2D8 = ([
        8 8 8 8 8 8 8 8 8 8
        8 8 8 8 8 8 7 7 7 8
        8 8 8 8 8 8 7 6 7 8
        8 8 8 8 8 8 7 7 7 8
        8 8 8 8 8 8 8 8 8 8])
    connectivity_2D8 = trues(3, 3)
    output_2D8 = hminima(image, connectivity_2D8, (3))
    @test output_2D8 == expected_connectivity_2D8
end

@testset "regional_maxima" begin
    image = ([
        0 0 0 0 0 0 0
        0 1 1 1 1 1 0
        0 1 2 2 2 1 0
        0 1 2 3 2 1 0
        0 1 2 2 2 1 0
        0 1 1 1 1 1 0
        0 0 0 0 0 0 0])

    expected_connectivity_2D8 = ([
        0 0 0 0 0 0 0
        0 0 0 0 0 0 0
        0 0 0 0 0 0 0
        0 0 0 1 0 0 0
        0 0 0 0 0 0 0
        0 0 0 0 0 0 0
        0 0 0 0 0 0 0])
    connectivity_2D8 = trues(3, 3)
    output_2D8 = regional_maxima(image, connectivity_2D8)
    @test output_2D8 == expected_connectivity_2D8
end

@testset "regional_minima" begin
    image = ([
        0 0 0 0 0 0 0
        0 3 3 3 3 3 0
        0 3 2 2 2 3 0
        0 3 2 1 2 3 0
        0 3 2 2 2 3 0
        0 3 3 3 3 3 0
        0 0 0 0 0 0 0])


    expected_connectivity_2D8 = ([
        1 1 1 1 1 1 1
        1 0 0 0 0 0 1
        1 0 0 0 0 0 1
        1 0 0 1 0 0 1
        1 0 0 0 0 0 1
        1 0 0 0 0 0 1
        1 1 1 1 1 1 1])
    connectivity_2D8 = trues(3, 3)
    output_2D8 = regional_minima(image, connectivity_2D8)
    @test output_2D8 == expected_connectivity_2D8
end
