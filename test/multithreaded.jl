# This doesn't run under "runtests.jl" but does under CI
using ImageMorphology
using Test

@testset "multithreaded" begin
    @test Threads.nthreads() > 1
    @testset "feature_transform" begin
        img = rand(100, 128) .> 0.9
        @test feature_transform(img; nthreads=Threads.nthreads()) == feature_transform(img; nthreads=1)
        # Since the threaded implementation handles two dimensions specially, we should check 0d and 1d
        img = reshape([true])
        @test feature_transform(img; nthreads=Threads.nthreads()) == feature_transform(img; nthreads=1) == reshape([CartesianIndex()])
        img = rand(100) .> 0.9
        @test feature_transform(img; nthreads=Threads.nthreads()) == feature_transform(img; nthreads=1)
    end
end
