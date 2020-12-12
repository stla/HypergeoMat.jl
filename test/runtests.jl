using HypergeoMat
using Test

@testset "Hypergeo2F1" begin
  @testset "some values" begin
    @test isapprox(
      hypergeomPQ(10, [1.0; 2.0], [3.0], [0.2; 0.5]),
      1.79412894456143
    )
    @test isapprox(
      hypergeomPQ(10, [1.0im; 2.0], [3.0im], [0.2; 0.5]),
      1.677558924 - 0.183004016im
    )
    @test isapprox(
      hypergeomPQ(10, [1.0; 2.0], [3.0], [0.2im; 0.5]),
      1.513810425 + 0.20576184im
    )
    @test isapprox(
      hypergeomPQ(10, [1.0; 2.0im], [3.0], [0.2im; 0.5]),
      0.7733140719 + 0.3092059749im
    )
  end
end
