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
  @testset "Gauss formula" begin
    a = 1.1
    b = 2.0im
    c = 9.0
    @test isapprox(
      hypergeomPQ(100, [a; b], [c], [1.0, 1.0, 1.0]),
      mvgamma(c, 3) * mvgamma(c-a-b, 3) / mvgamma(c-a, 3) / mvgamma(c-b, 3)
    )
  end
end
