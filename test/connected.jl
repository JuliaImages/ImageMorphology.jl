@testset "Label components" begin
    A = [true  true  false true;
         true  false true  true]
    lbltarget = [1 1 0 2;
                 1 0 2 2]
    lbltarget1 = [1 2 0 4;
                  1 0 3 4]
    @test label_components(A) == lbltarget
    @test label_components(A, [1]) == lbltarget1
    connectivity = [false true  false;
                    true  false true;
                    false true  false]
    @test label_components(A, connectivity) == lbltarget
    connectivity = trues(3,3)
    lbltarget2 = [1 1 0 1;
                  1 0 1 1]
    @test label_components(A, connectivity) == lbltarget2
end
