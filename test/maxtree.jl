@testset "connectivity neighborhoods" begin
    @test @inferred(ImageMorphology.neighbor_cartesian_offsets(CartesianIndex{1}, 1)) ==
          ImageMorphology.neighbor_cartesian_offsets(CartesianIndex{1}, 2) ==
          [CartesianIndex{1}(-1), CartesianIndex{1}(1)]
    conn2d_4 = ImageMorphology.neighbor_cartesian_offsets(CartesianIndex{2}, 1)
    @test conn2d_4 ==
        [CartesianIndex{2}(0, -1), CartesianIndex{2}(-1, 0),
         CartesianIndex{2}(1, 0), CartesianIndex{2}(0, 1)]
    @test @inferred(ImageMorphology.neighbor_cartesian_offsets(CartesianIndex{2}, 2)) ==
          [CartesianIndex{2}(-1, -1), CartesianIndex{2}(0, -1), CartesianIndex{2}(1, -1),
           CartesianIndex{2}(-1, 0), CartesianIndex{2}(1, 0),
           CartesianIndex{2}(-1, 1), CartesianIndex{2}(0, 1), CartesianIndex{2}(1, 1)]
    @test ImageMorphology.linear_offsets(conn2d_4, fill(0, (1, 1))) == [-1, -1, 1, 1]
    @test ImageMorphology.linear_offsets(conn2d_4, fill(0, (2, 2))) == [-2, -1, 1, 2]
    @test @inferred(ImageMorphology.linear_offsets(conn2d_4, fill(0, (3, 3)))) == [-3, -1, 1, 3]
    @test @inferred(ImageMorphology.linear_offsets(conn2d_4, fill(0, (7, 9)))) == [-7, -1, 1, 7]
end

const CI2 = CartesianIndex{2}

bbox2d(x1, y1, x2, y2) = (CI2(x1, y1), CI2(x2, y2))

@testset "MaxTree construction" begin
    A = [15 13 16;
         11 12 10;
         16 11 14]
    mtree = MaxTree(A, connectivity=1, rev=false) # test kwargs defaults
    @test mtree isa MaxTree{2}
    @test mtree == MaxTree(A) # test kwargs defaults
    @test isequal(mtree, MaxTree(A))
    @test ndims(mtree) == 2
    @test size(mtree) == size(A)
    @test axes(mtree) == axes(A)
    @test length(mtree) == length(A)
    @test !mtree.rev
    @test mtree.parentindices == [4 5 4;
                                  8 2 8;
                                  2 2 2]
    @test mtree.traverse == [8, 2, 6, 5, 4, 9, 1, 3, 7]
    @test areas(mtree) == [1 3 1; 8 4 9; 1 1 1]
    @test diameters(mtree) ==[1 3 1; 3 3 3; 1 1 1]
    @test boundingboxes(mtree) ==
        [bbox2d(1, 1, 1, 1) bbox2d(1, 1, 1, 3) bbox2d(1, 3, 1, 3);
         bbox2d(1, 1, 3, 3) bbox2d(1, 1, 2, 3) bbox2d(1, 1, 3, 3);
         bbox2d(3, 1, 3, 1) bbox2d(3, 2, 3, 2) bbox2d(3, 3, 3, 3)]

    # test 8-neighborhood
    mtree2 = MaxTree(A, connectivity=2)
    @test !mtree2.rev
    @test mtree2 != mtree
    @test mtree2.parentindices == [4 5 4;
                                   8 2 8;
                                   5 2 5]
    @test mtree2.traverse == [8, 2, 6, 5, 4, 9, 1, 3, 7]
    @test areas(mtree2) == [1 3 1; 8 6 9; 1 1 1]
    @test diameters(mtree2) == [1 3 1; 3 3 3; 1 1 1]
    @test boundingboxes(mtree2) ==
        [bbox2d(1, 1, 1, 1) bbox2d(1, 1, 1, 3) bbox2d(1, 3, 1, 3);
         bbox2d(1, 1, 3, 3) bbox2d(1, 1, 3, 3) bbox2d(1, 1, 3, 3);
         bbox2d(3, 1, 3, 1) bbox2d(3, 2, 3, 2) bbox2d(3, 3, 3, 3)]

    # test reverse (brightest to darkest) MaxTree
    mtree_rev = MaxTree(A, rev=true)
    @test mtree_rev.rev
    A = [15 13 16;
         11 12 10;
         16 11 14]
    @test mtree_rev.parentindices == [3 9 3;
                                      5 4 5;
                                      3 5 1]
    @test mtree_rev.traverse == [3, 7, 1, 9, 4, 5, 2, 6, 8]
    @test areas(mtree_rev) == [7 5 1; 1 4 1; 9 1 6]
    @test diameters(mtree_rev) == [3 3 1; 1 3 1; 3 1 3]
    @test boundingboxes(mtree_rev) ==
        [bbox2d(1, 1, 3, 3) bbox2d(1, 1, 3, 3) bbox2d(1, 3, 1, 3);
         bbox2d(2, 1, 2, 1) bbox2d(2, 1, 3, 3) bbox2d(2, 3, 2, 3);
         bbox2d(1, 1, 3, 3) bbox2d(3, 2, 3, 2) bbox2d(1, 1, 3, 3)]

    # test reverse 8-neighborhood MaxTree
    mtree2_rev = MaxTree(A, rev=true, connectivity=2)
    @test mtree2_rev.rev
    @test mtree2_rev.parentindices == [3 9 3;
                                       5 4 2;
                                       3 2 1]
    @test mtree2_rev.traverse == [3, 7, 1, 9, 4, 5, 2, 6, 8]
    @test areas(mtree2_rev) == [7 5 1; 3 4 1; 9 1 6]
    @test diameters(mtree2_rev) == [3 3 1; 3 3 1; 3 1 3]
    @test boundingboxes(mtree2_rev) ==
        [bbox2d(1, 1, 3, 3) bbox2d(1, 1, 3, 3) bbox2d(1, 3, 1, 3)
         bbox2d(2, 1, 3, 3) bbox2d(2, 1, 3, 3) bbox2d(2, 3, 2, 3)
         bbox2d(1, 1, 3, 3) bbox2d(3, 2, 3, 2) bbox2d(1, 1, 3, 3)]

    # degenerated cases
    A00 = fill(0, (0, 0))
    mtree00 = MaxTree(A00, connectivity=1, rev=false)
    @test mtree00 isa MaxTree{2}
    @test areas(mtree00) == Matrix{Int}(undef, 0, 0)

    A11 = fill(0, (1, 1))
    mtree11 = MaxTree(A11, connectivity=2, rev=false)
    @test mtree11 isa MaxTree{2}
    @test areas(mtree11) == fill(1, (1, 1))

    A101 = fill(0, (1, 0, 1))
    mtree101 = MaxTree(A101, connectivity=2, rev=false)
    @test mtree101 isa MaxTree{3}
    @test diameters(mtree101) == Array{Int}(undef, 1, 0, 1)
end

@testset "local_minima/maxima()" begin
    A = [15 13 16;
         11 12 10;
         16 11 14]
    @test local_maxima(A) == [3 0 1;
                              0 0 0;
                              2 0 4]
    B = similar(A, Float64)
    @test local_maxima!(B, A) === B
    @test B == float(local_maxima(A))

    @test_throws DimensionMismatch local_maxima!(fill(0, (2, 2)), A)

    @test local_maxima(A, connectivity=2) == [3 0 1; 0 0 0; 2 0 4]
    atree = MaxTree(A)
    atree_rev = MaxTree(A, rev=true)
    @test local_maxima(A, maxtree=atree, connectivity=2) == local_maxima(A)
    @test_throws ArgumentError local_maxima(A, maxtree=atree_rev)

    # test local_minima
    @test local_minima(A) == [0 0 0;
                              3 0 1;
                              0 2 0]
    @test local_minima(A, connectivity=2) == [0 0 0;
                                              0 0 1;
                                              0 0 0]
    @test local_minima!(B, A) === B
    @test B == float(local_minima(A))
    @test_throws DimensionMismatch local_minima!(fill(0, (4, 2)), A)

    @test local_minima(A, maxtree=atree_rev, connectivity=2) == local_minima(A)
    @test_throws ArgumentError local_minima(A, maxtree=atree)
end

@testset "area_opening/closing()" begin
    A = [3 1 3 1 4 4;
         1 3 1 1 3 3;
         1 1 1 1 1 1;
         1 4 4 1 1 2;
         1 4 1 1 5 5;
         1 1 1 2 5 2]
    B = similar(A)
    @test B === area_opening!(B, A)
    @test_throws DimensionMismatch area_opening!(similar(A, (7, 6)), A)
    atree = MaxTree(A)
    atree_rev = MaxTree(A, rev=true)
    atree_small = MaxTree(A[1:5, 2:6])
    @test B == area_opening(A, maxtree=atree)
    @test_throws ArgumentError area_opening(A, maxtree=atree_rev)
    @test_throws DimensionMismatch area_opening(A, maxtree=atree_small)
    @test area_opening(A, min_area=3) ==
        [1 1 1 1 3 3;
         1 1 1 1 3 3;
         1 1 1 1 1 1;
         1 4 4 1 1 2;
         1 4 1 1 5 5;
         1 1 1 2 5 2]
    @test area_opening(A, min_area=1) == A
    @test area_opening(A, min_area=5) ==
        [1 1 1 1 1 1;
         1 1 1 1 1 1;
         1 1 1 1 1 1;
         1 1 1 1 1 2;
         1 1 1 1 2 2;
         1 1 1 2 2 2]
    @test area_opening(A, min_area=45) == fill(0, (6, 6)) # image less than min_area

    mA = .-A
    mB = similar(mA)
    @test mB === area_closing!(mB, mA)
    @test_throws DimensionMismatch area_closing!(similar(A, (7, 6)), A)
    matree = MaxTree(mA)
    matree_rev = MaxTree(mA, rev=true)
    @test mB == area_closing(mA, maxtree=atree_rev)
    @test_throws ArgumentError area_closing(mA, maxtree=matree)
    @test area_closing(mA, min_area=3) == .-area_opening(A, min_area=3)
end

@testset "diameter_opening/closing()" begin
    A = [3 1 3 1 3 4;
         1 3 1 1 4 3;
         1 1 4 1 1 3;
         1 4 4 1 1 2;
         1 4 1 1 1 2;
         1 1 1 2 5 1]
    B = similar(A)
    @test B === diameter_opening!(B, A)
    @test_throws DimensionMismatch diameter_opening!(similar(A, (7, 6)), A)
    atree = MaxTree(A)
    atree_rev = MaxTree(A, rev=true)
    atree_small = MaxTree(A[1:5, 2:6])
    @test B == diameter_opening(A, maxtree=atree)
    @test_throws ArgumentError diameter_opening(A, maxtree=atree_rev)
    @test_throws DimensionMismatch diameter_opening(A, maxtree=atree_small)
    @test diameter_opening(A, min_diameter=3) ==
        [1 1 1 1 3 3;
         1 1 1 1 3 3;
         1 1 4 1 1 3;
         1 4 4 1 1 2;
         1 4 1 1 1 2;
         1 1 1 1 1 1]
    @test diameter_opening(A, min_diameter=1) == A
    @test diameter_opening(A, min_diameter=3, connectivity=2) ==
        [3 1 3 1 3 3;
         1 3 1 1 3 3;
         1 1 4 1 1 3;
         1 4 4 1 1 2;
         1 4 1 1 1 2;
         1 1 1 2 2 1]
     @test diameter_opening(A, min_diameter=5, connectivity=2) ==
        [3 1 3 1 2 2;
         1 3 1 1 2 2;
         1 1 3 1 1 2;
         1 3 3 1 1 2;
         1 3 1 1 1 2;
         1 1 1 2 2 1]
    @test diameter_opening(A, min_diameter=45) == fill(0, (6, 6)) # image less than min_diameter

    mA = .-A
    mB = similar(mA)
    @test mB === diameter_closing!(mB, mA)
    @test_throws DimensionMismatch diameter_closing!(similar(A, (7, 6)), A)
    matree = MaxTree(mA)
    matree_rev = MaxTree(mA, rev=true)
    @test mB == diameter_closing(mA, maxtree=atree_rev)
    @test_throws ArgumentError diameter_closing(mA, maxtree=matree)
    @test diameter_closing(mA, min_diameter=3) == .-diameter_opening(A, min_diameter=3)
end
