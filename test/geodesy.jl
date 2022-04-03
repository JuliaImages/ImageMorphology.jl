@testset "underbuild" begin
    marker = ([0 0 0 0 0;
              0 9 0 0 0;
              0 0 0 0 0;
              0 0 0 3 0;
              0 0 0 0 0;
              0 9 0 0 0])

    mask = ([0 0 0 0 0;
            0 8 7 1 0;
            0 9 0 4 0;
            0 0 0 4 0;
            0 0 6 5 6;
            0 0 9 8 9])

    expected_connectivity_2D8 = ([0 0 0 0 0;
                                 0 8 7 1 0;
                                 0 8 0 4 0;
                                 0 0 0 4 0;
                                 0 0 4 4 4;
                                 0 0 4 4 4])  
    connectivity_2D8 = trues(3,3)
    output_2D8 = underbuild(marker,mask,connectivity_2D8)
    @test output_2D8 == expected_connectivity_2D8

    expected_connectivity_2D4 = ([0 0 0 0 0;
            0 8 7 1 0;
            0 8 0 3 0;
            0 0 0 3 0;
            0 0 3 3 3;
            0 0 3 3 3])
    connectivity_2D4 = [false true  false;
                    true  false true;
                    false true  false]
    output_2D4 = underbuild(marker,mask,connectivity_2D4)
    @test output_2D4 == expected_connectivity_2D4
end

@testset "overbuild" begin
    marker = ([9 9 9 9 9;
    9 8 9 9 9;
    9 9 9 9 9;
    9 9 9 3 9;
    9 9 9 9 9;
    9 8 9 9 9])

    mask = ([9 9 9 9 9;
    9 2 3 8 9;
    9 1 9 6 9;
    9 9 9 6 9;
    9 9 4 5 4;
    9 9 1 2 1])

    expected_connectivity_2D8 = ([9 9 9 9 9;
    9 6 6 8 9;
    9 6 9 6 9;
    9 9 9 6 9;
    9 9 6 6 6;
    9 9 6 6 6])  
    connectivity_2D8 = trues(3,3)
    output_2D8 = overbuild(marker,mask,connectivity_2D8)
    @test output_2D8 == expected_connectivity_2D8

    expected_connectivity_2D4 = ([9 9 9 9 9;
    9 8 8 8 9;
    9 8 9 6 9;
    9 9 9 6 9;
    9 9 6 6 6;
    9 9 6 6 6])
    connectivity_2D4 = [false true  false;
                    true  false true;
                    false true  false]
    output_2D4 = overbuild(marker,mask,connectivity_2D4)
    @test output_2D4 == expected_connectivity_2D4
end

@testset "hmaxima" begin
    image = ([
    3 3 3 3 3 3 3 3 3 3;
    3 4 4 4 3 3 4 4 4 3;
    3 4 5 4 3 3 4 9 4 3;
    3 4 4 4 3 3 4 4 4 3;
    3 3 3 3 3 3 3 3 3 3])

    expected_connectivity_2D8 =([
    3 3 3 3 3 3 3 3 3 3;
    3 3 3 3 3 3 4 4 4 3;
    3 3 3 3 3 3 4 6 4 3;
    3 3 3 3 3 3 4 4 4 3;
    3 3 3 3 3 3 3 3 3 3])
    connectivity_2D8 = trues(3,3)
    output_2D8 = hmaxima(image,connectivity_2D8,(3))
    @test output_2D8 == expected_connectivity_2D8
end
