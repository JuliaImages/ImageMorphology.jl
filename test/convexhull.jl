@testset "Convex Hull" begin
    A = zeros(50, 30)
    A = convert(Array{Bool}, A)
    A[25, 1] = 1
    A[1, 10] = 1
    A[10, 10] = 1
    A[10, 30] = 1
    A[40, 30] = 1
    A[40, 10] = 1
    A[50, 10] = 1
    B = @inferred convexhull(A)
    C = CartesianIndex{}[]
    push!(C, CartesianIndex{}(25, 1))
    push!(C, CartesianIndex{}(1, 10))
    push!(C, CartesianIndex{}(10, 30))
    push!(C, CartesianIndex{}(40, 30))
    push!(C, CartesianIndex{}(50, 10))
    @test typeof(B) == Array{CartesianIndex{2},1}
    @test sort(B) == sort(C)

    A = [
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
        0.0,
        1.0,
        1.0,
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
        1.0,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
    ]
    A = reshape(A, 5, 5)
    A = convert(Array{Bool}, A)
    B = B = @inferred convexhull(A)
    C = CartesianIndex{}[]
    push!(C, CartesianIndex{}(1, 3))
    push!(C, CartesianIndex{}(3, 1))
    push!(C, CartesianIndex{}(3, 5))
    push!(C, CartesianIndex{}(5, 3))
    @test typeof(B) == Array{CartesianIndex{2},1}
    @test sort(B) == sort(C)
end
