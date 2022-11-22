########################
# Analytic Image Tests #
########################
using Test: @testset, @test
using Apolo.Images
using Apolo.Images.AnalyticImages

const INTERVAL_START = LinRange(-10.0, 10.0, 20)
const INTERVAL_POS = LinRange(1.0, 10.0, 20)
const INTERVAL_OFFSET = LinRange(-10.0, -1.0, 20)
const TOLERANCE = 1e-3

@testset "Images.AnalyticImages 3D" begin

    analytic_intensity(x, y, z) = sin(y) * cos(x) * tan(z)
    start_img = (rand(INTERVAL_START), rand(INTERVAL_START), rand(INTERVAL_START))
    length_img = (rand(INTERVAL_POS), rand(INTERVAL_POS), rand(INTERVAL_POS))
    offset_img = (0.2, 0.1, 0.3)
    finish_img = start_img .+ length_img

    a_img = AnalyticImage(analytic_intensity, start_img, finish_img)

    @test start(a_img) == start_img
    @test finish(a_img) == finish_img
    @test length(a_img) == finish_img .- start_img
    @test dimension(a_img) == length(start_img) == length(finish_img)
    @test intensity(a_img) == analytic_intensity
    @test grid(a_img) == nothing

    @test a_img((@. start_img + length_img / 100 - offset_img)..., offset=offset_img) ≈
          analytic_intensity((@. start_img + length_img / 100)...) atol = TOLERANCE
    @test a_img((@. finish_img - length_img / 100 - offset_img)..., offset=offset_img) ≈
          analytic_intensity((@. finish_img - length_img / 100)...) atol = TOLERANCE
    vec_points = [start_img, finish_img]
    @test a_img(vec_points) ≈
          [analytic_intensity(start_img...), analytic_intensity(finish_img...)] atol = TOLERANCE

end
