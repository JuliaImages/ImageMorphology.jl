@testset "low_leveling" begin
    ref = ([
        2 2 2;
        5 5 5;
        1 1 1])

    marker = ([
        1 1 1;
        3 3 3;
        4 4 4])

    expected_connectivity_2D8 = ([
        2 2 2;
        4 4 4;
        1 1 1])  
    connectivity_2D8 = trues(3,3)
    output_2D8 = low_leveling(ref,marker,connectivity_2D8)
    @test output_2D8 == expected_connectivity_2D8

    ref = ([
        5 0 0 0 0 0 5;
        0 5 0 0 0 5 0;
        0 0 5 5 5 0 0;
        0 0 5 5 5 0 0;
        0 0 5 5 5 0 0;
        0 5 0 0 0 5 0;
        5 0 0 0 0 0 5])

    marker = ([
        0 0 0 0 0 0 0;
        0 0 0 0 0 0 0;
        0 0 6 6 6 0 0;
        0 0 6 6 6 0 0;
        0 0 6 6 6 0 0;
        0 0 0 0 0 0 0;
        0 0 0 0 0 0 0])

    expected_connectivity_2D8 = ([
        5 0 0 0 0 0 5;
        0 5 0 0 0 5 0;
        0 0 5 5 5 0 0;
        0 0 5 5 5 0 0;
        0 0 5 5 5 0 0;
        0 5 0 0 0 5 0;
        5 0 0 0 0 0 5])  
    connectivity_2D8 = trues(3,3)
    output_2D8 = low_leveling(ref,marker,connectivity_2D8)
    @test output_2D8 == expected_connectivity_2D8

end

@testset "high_leveling" begin
    ref = ([
        2 2 2;
        5 5 5;
        1 1 1])

    marker = ([
        1 1 1;
        3 3 3;
        4 4 4])

    expected_connectivity_2D8 = ([
        2 2 2; 
        5 5 5; 
        3 3 3])  
    connectivity_2D8 = trues(3,3)
    output_2D8 = high_leveling(ref,marker,connectivity_2D8)
    @test output_2D8 == expected_connectivity_2D8

    ref = ([
        5 0 0 0 0 0 5;
        0 5 0 0 0 5 0;
        0 0 5 5 5 0 0;
        0 0 5 5 5 0 0;
        0 0 5 5 5 0 0;
        0 5 0 0 0 5 0;
        5 0 0 0 0 0 5])

    marker = ([
        0 0 0 0 0 0 0;
        0 0 0 0 0 0 0;
        0 0 4 4 4 0 0;
        0 0 4 4 4 0 0;
        0 0 4 4 4 0 0;
        0 0 0 0 0 0 0;
        0 0 0 0 0 0 0])

    expected_connectivity_2D8 = ([
        5 0 0 0 0 0 5;
        0 5 0 0 0 5 0;
        0 0 5 5 5 0 0;
        0 0 5 5 5 0 0;
        0 0 5 5 5 0 0;
        0 5 0 0 0 5 0;
        5 0 0 0 0 0 5])  
    connectivity_2D8 = trues(3,3)
    output_2D8 = high_leveling(ref,marker,connectivity_2D8)
    @test output_2D8 == expected_connectivity_2D8

end

@testset "leveling" begin
    
    ref = ([1 1 2 1 1 4 4 4 4 4 0 0 0 0 0])
    marker = ([0 0 0 0 0 0 6 6 6 6 6 0 0 2 0])
    expected = ([1 1 1 1 1 4 4 4 4 4 0 0 0 0 0])
    connectivity = trues(1,3)
    output = leveling(ref,marker,connectivity)
    @test output == expected
    
    ref= ([1 1 1 1 5 5 5 5])
    marker=([4 4 4 4 2 2 2 2])
    expected=([2 2 2 2 4 4 4 4])
    connectivity = trues(1,3)
    output = leveling(ref,marker,connectivity)
    @test output == expected

end