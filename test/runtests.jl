using HypergeoMat
using Test
using LinearAlgebra

@testset "Hypergeo0F0 is exponential of trace" begin
  X1 = [0.3 0.2 0.1; 0.2 0.3 0.2; 0.1 0.2 0.3]
  X2 = [0.3im 0.2 0.1; 0.2 0.3im 0.2; 0.1 0.2 0.3im]
  @test isapprox(
    hypergeomPQ(10, Float64[], Float64[], X1),
    exp(tr(X1))
  )
  @test isapprox(
    hypergeomPQ(10, Float64[], Float64[], X2),
    exp(tr(X2))
  )
end

@testset "Hypergeo1F0 is det(I-X)^(-a)" begin
  X1 = [0.03 0.02 0.01; 0.02 0.03 0.02; 0.01 0.02 0.03]
  X2 = [0.02 0.02 0.01; 0.02 0.02 0.02; 0.01 0.02 0.02]
  X3 = [0.03im 0.02 0.01; 0.02 0.03im 0.02; 0.01 0.02 0.03im]
  @test isapprox(
    hypergeomPQ(15, [3.0], Float64[], X1),
    det(I - X1)^(-3)
  )
  @test isapprox(
    hypergeomPQ(15, [4.0im], Float64[], X2),
    det(I - X2)^(-4.0im)
  )
  @test isapprox(
    hypergeomPQ(15, [3.0], Float64[], X3),
    det(I - X3)^(-3)
  )
end

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
