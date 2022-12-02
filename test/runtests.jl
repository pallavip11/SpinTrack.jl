using SpinTrack
using Test
using BenchmarkTools

@testset "SpinTrack.jl" begin
    # Write your tests here.
end

@testset "Performance" begin
    @benchmark get_solution(u0_long, symmetric_hybrid_ring())

    @test (1.0 > @elapsed get_solution(u0_long, symmetric_hybrid_ring()))
end
