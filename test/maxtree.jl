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
    @test mtree.parents == [4 5 4;
                            8 2 8;
                            2 2 2]
    @test mtree.traverse == [8, 2, 6, 5, 4, 9, 1, 3, 7]
    @test areas(mtree) == [1, 8, 1, 3, 4, 1, 1, 9, 1]
    @test diameters(mtree) == [1, 3, 1, 3, 3, 1, 1, 3, 1]
    @test boundingboxes(mtree) == [1 1 3 1 1 3 1 1 3;
                                   1 1 1 1 1 2 3 1 3;
                                   1 3 3 1 2 3 1 3 3;
                                   1 3 1 3 3 2 3 3 3]
    # test 8-neighborhood
    mtree2 = MaxTree(A, connectivity=2)
    @test !mtree2.rev
    @test mtree2 != mtree
    @test mtree2.parents == [4 5 4;
                             8 2 8;
                             5 2 5]
    @test mtree2.traverse == [8, 2, 6, 5, 4, 9, 1, 3, 7]
    @test areas(mtree2) == [1, 8, 1, 3, 6, 1, 1, 9, 1]
    @test diameters(mtree2) == [1, 3, 1, 3, 3, 1, 1, 3, 1]
    @test boundingboxes(mtree2) == [1 1 3 1 1 3 1 1 3;
                                    1 1 1 1 1 2 3 1 3;
                                    1 3 3 1 3 3 1 3 3;
                                    1 3 1 3 3 2 3 3 3]

    # test reverse (brightest to darkest) MaxTree
    mtree_rev = MaxTree(A, rev=true)
    @test mtree_rev.rev
    A = [15 13 16;
         11 12 10;
         16 11 14]
    @test mtree_rev.parents == [3 9 3;
                                5 4 5;
                                3 5 1]
    @test mtree_rev.traverse == [3, 7, 1, 9, 4, 5, 2, 6, 8]
    @test areas(mtree_rev) == [7, 1, 9, 5, 4, 1, 1, 1, 6]
    @test diameters(mtree_rev) == [3, 1, 3, 3, 3, 1, 1, 1, 3]
    @test boundingboxes(mtree_rev) == [1 2 1 1 2 3 1 2 1;
                                       1 1 1 1 1 2 3 3 1;
                                       3 2 3 3 3 3 1 2 3;
                                       3 1 3 3 3 2 3 3 3]

    # test reverse 8-neighborhood MaxTree
    mtree2_rev = MaxTree(A, rev=true, connectivity=2)
    @test mtree2_rev.rev
    @test mtree2_rev.parents == [3 9 3;
                                 5 4 2;
                                 3 2 1]
    @test mtree2_rev.traverse == [3, 7, 1, 9, 4, 5, 2, 6, 8]
    @test areas(mtree2_rev) == [7, 3, 9, 5, 4, 1, 1, 1, 6]
    @test diameters(mtree2_rev) == [3, 3, 3, 3, 3, 1, 1, 1, 3]
    @test boundingboxes(mtree2_rev) == [1 2 1 1 2 3 1 2 1;
                                        1 1 1 1 1 2 3 3 1;
                                        3 3 3 3 3 3 1 2 3;
                                        3 3 3 3 3 2 3 3 3]
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
