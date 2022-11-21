#######################
# Ferrite Image Tests #
#######################
using Test: @testset, @test
using Apolo.Geometry.FerriteGrids: _convert_to_ferrite_nomenclature
using Apolo.Images
using Apolo.Images.FerriteImages
using Apolo.Images: _total_num_pixels, _index_is_inbounds, _eval_intensity,
    _interpolate, _extrapolate

const INTERVAL_START = LinRange(-10.0, 10.0, 20)
const INTERVAL_POS = LinRange(1.0, 10.0, 20)
const INTERVAL_OFFSET = LinRange(-10.0, -1.0, 20)
const TOLERANCE = 1e-3

@testset "Images.FerriteImage 2D" begin

    start_img = (rand(INTERVAL_START), rand(INTERVAL_START))
    spacing_img = (rand(INTERVAL_POS), rand(INTERVAL_POS))
    num_pixels_img = (4, 3)
    start_img_grid = start_img .+ spacing_img ./ 2
    image_dimension = length(spacing_img)
    finish_img = start_img .+ num_pixels_img .* spacing_img
    finish_img_grid = finish_img .- spacing_img ./ 2
    length_img = finish_img .- start_img
    @test finish(start_img, num_pixels_img, spacing_img) == finish_img
    @test length(start_img, finish_img) == length_img

    # image intensity
    intensity_vec = collect(1.0:Float64(prod(num_pixels_img)))
    intensity_array = reshape(intensity_vec, num_pixels_img)

    # generate default image
    f_img_default = FerriteImage(intensity_array)
    @test start(f_img_default) == (0, 0)
    @test spacing(f_img_default) == (1, 1)
    @test finish(f_img_default) == size(intensity(f_img_default))
    @test length(f_img_default) == finish(f_img_default)

    # generate image and internaly cash a grid
    f_img = FerriteImage(intensity_array, spacing_img, start_img)

    fgrid_img = grid(f_img)
    @test grid(f_img).grid.cells == fgrid_img.grid.cells

    # generate a ferrite image with a given grid
    # f_img = FerriteImage(intensity_array, fgrid_img, spacing_img, start_img)

    # test the converter intensity function
    fintensity_to_test, cart_indexes_to_test = _convert_to_ferrite_nomenclature(
        intensity_array,
        fgrid_img
    )
    fintensity_hand = [1, 2, 6, 5, 3, 7, 4, 8, 10, 9, 11, 12] # this must be the one dim pixel index
    @test fintensity_to_test ≈ fintensity_hand atol = TOLERANCE
    cartesian_indexes_hand = cartesian_index.(fintensity_hand, Ref(fgrid_img))
    @test cart_indexes_to_test == cartesian_indexes_hand

    # test getter functions
    @test start(f_img) == start_img
    @test start_grid(f_img) == start(grid(f_img)) == start_img_grid
    @test finish(f_img) == finish_img
    @test finish_grid(f_img) == finish(grid(f_img)) == finish_img_grid
    @test collect(length(f_img)) ≈ collect(length_img) atol = TOLERANCE
    coords = coordinates(f_img)
    @test coords == LinRange.(start_img_grid, finish_img_grid, num_pixels_img)
    @test spacing_img[1] ≈ (coords[1][2] - coords[1][1]) atol = TOLERANCE
    @test spacing_img[2] ≈ (coords[2][2] - coords[2][1]) atol = TOLERANCE

    @test extrema(f_img) == [(start_img[1], finish_img[1]), (start_img[2], finish_img[2])]
    @test intensity(f_img) == intensity_array
    @test intensity_type(f_img) == intensity_type(f_img.intensity)
    @test dimension(f_img) == image_dimension
    @test num_pixels(f_img) == num_pixels_img
    @test _total_num_pixels(f_img) == prod(num_pixels_img)
    @test spacing(f_img) == spacing_img

    # test checker functions
    # ∈ and ∉
    @test start_img ∈ f_img
    @test finish_img ∈ f_img
    @test (start_img .+ (finish_img .- start_img) ./ 2) ∈ f_img
    @test (start_img .+ 9 .* finish_img) ∉ f_img

    # indexes
    # is inside function
    @test !(_index_is_inbounds(-1, f_img))
    @test !(_index_is_inbounds(_total_num_pixels(f_img) + 1, f_img))
    # normal case
    @test cartesian_index(10, f_img) == CartesianIndex(2, 3)
    # edge cases
    @test cartesian_index(2, f_img) == CartesianIndex(2, 1)
    @test cartesian_index(12, f_img) == CartesianIndex(4, 3)

    # test interpolation functions
    @test _interpolate([start_img_grid], f_img) ≈ [intensity_array[1, 1]] atol = TOLERANCE
    @test [f_img(start_img_grid...)] ≈ [intensity_array[1, 1]] atol = TOLERANCE

    # final point the image grid
    @test _extrapolate([finish_img], f_img) ≈ [intensity_array[end, end]] atol = TOLERANCE
    @test [intensity_array[end, end]] ≈ [f_img(finish_img...)] atol = TOLERANCE

    # test point at the final cell middle point
    final_mid_cell_point = finish_img_grid .- spacing_img ./ 2
    hand_intensity = sum(intensity_array[end-1:end, end-1:end]) / 4
    @test _interpolate([final_mid_cell_point], f_img) ≈ [hand_intensity] atol = TOLERANCE
    @test [hand_intensity] ≈ [f_img(final_mid_cell_point...)] atol = TOLERANCE

    # test functor with a vector of points
    vec_points = [start_img_grid, finish_img, final_mid_cell_point]
    vec_intensities_hand = [intensity_array[1, 1], intensity_array[end, end], hand_intensity]
    @test _eval_intensity(vec_points, f_img) ≈ f_img(vec_points) atol = TOLERANCE
    @test f_img(vec_points) ≈ vec_intensities_hand atol = TOLERANCE

end

@testset "Images.FerriteImage 3D" begin

    start_img = (rand(INTERVAL_START), rand(INTERVAL_START), rand(INTERVAL_START))
    spacing_img = (rand(INTERVAL_POS), rand(INTERVAL_POS), rand(INTERVAL_POS))
    num_pixels_img = (2, 4, 3)
    start_img_grid = start_img .+ spacing_img ./ 2
    image_dimension = length(spacing_img)
    finish_img = start_img .+ num_pixels_img .* spacing_img
    finish_img_grid = finish_img .- spacing_img ./ 2
    length_img = finish_img .- start_img

    # image intensity
    intensity_vec = collect(1.0:Float64(prod(num_pixels_img)))
    intensity_array = reshape(intensity_vec, num_pixels_img)

    # create the image ferrite grid
    fgrid_img = create_ferrite_img_fgrid(start_img, spacing_img, length_img, num_pixels_img)

    # test the converter intensity function
    fintensity_to_test, cart_indexes_to_test = _convert_to_ferrite_nomenclature(intensity_array, fgrid_img)

    fintensity_hand = [
        1, 2, 4, 3, 9, 10, 12, 11, 6, 5, 14, 13, 8,
        7, 16, 15, 17, 18, 20, 19, 22, 21, 24, 23
    ] # this must be the one dim pixel index

    @test fintensity_to_test == fintensity_hand
    cartesian_indexes_hand = cartesian_index.(fintensity_hand, Ref(fgrid_img))
    @test fintensity_to_test == fintensity_hand
    @test cart_indexes_to_test == cartesian_indexes_hand

    # generate image and cash a grid
    f_img = FerriteImage(intensity_array, spacing_img, start_img)
    fgrid_img = grid(f_img)
    @test grid(f_img).grid.cells == fgrid_img.grid.cells

    # generate a ferrite image with a given grid
    # f_img = FerriteImage(intensity_array, fgrid_img, spacing_img, start_img)

    # test getter functions
    @test start(f_img) == start_img
    @test start_grid(f_img) == start(grid(f_img)) == start_img_grid
    @test finish(f_img) == finish_img
    @test finish_grid(f_img) == finish(grid(f_img)) == finish_img_grid
    @test collect(length(f_img)) ≈ collect(length_img) atol = TOLERANCE
    coords = coordinates(f_img)
    @test coords == LinRange.(start_img_grid, finish_img_grid, num_pixels_img)
    @test spacing_img[1] ≈ (coords[1][2] - coords[1][1]) atol = TOLERANCE
    @test spacing_img[2] ≈ (coords[2][2] - coords[2][1]) atol = TOLERANCE
    @test spacing_img[3] ≈ (coords[3][2] - coords[3][1]) atol = TOLERANCE

    @test extrema(f_img) == [
        (start_img[1], finish_img[1]),
        (start_img[2], finish_img[2]),
        (start_img[3], finish_img[3]),
    ]
    @test intensity(f_img) == intensity_array
    @test intensity_type(f_img) == intensity_type(f_img.intensity)
    @test dimension(f_img) == image_dimension
    @test num_pixels(f_img) == num_pixels_img
    @test _total_num_pixels(f_img) == prod(num_pixels_img)
    @test spacing(f_img) == spacing_img

    # test checker functions
    # ∈ and ∉
    @test start_img ∈ f_img
    @test finish_img ∈ f_img
    @test (start_img .+ (finish_img .- start_img) ./ 2) ∈ f_img
    @test (start_img .+ 900 .* finish_img) ∉ f_img

    # indexes
    # is inside function
    @test !(_index_is_inbounds(-1, f_img))
    @test !(_index_is_inbounds(_total_num_pixels(f_img) + 1, f_img))
    # normal case
    @test cartesian_index(10, f_img) == CartesianIndex(2, 1, 2)
    # edge cases
    @test cartesian_index(2, f_img) == CartesianIndex(2, 1, 1)
    @test cartesian_index(16, f_img) == CartesianIndex(2, 4, 2)

    # test interpolation functions
    @test _interpolate([start_img_grid], f_img) ≈ [intensity_array[1, 1, 1]] atol = TOLERANCE
    @test [f_img(start_img_grid...)] ≈ [intensity_array[1, 1, 1]] atol = TOLERANCE

    # final point the image grid
    @test _extrapolate([finish_img .+ (1.0, 1.0, 1.0)], f_img) ≈ [intensity_array[end, end, end]] atol = TOLERANCE
    @test [intensity_array[end, end, end]] ≈ [f_img(finish_img...)] atol = TOLERANCE

    # test point at the final cell middle point
    final_mid_cell_point = finish_img_grid .- spacing_img ./ 2
    hand_intensity = sum(intensity_array[end-1:end, end-1:end, end-1:end]) / 8
    @test _interpolate([final_mid_cell_point], f_img) ≈ [hand_intensity] atol = TOLERANCE
    @test [hand_intensity] ≈ [f_img(final_mid_cell_point...)] atol = TOLERANCE

    # test functor with a vector of points
    offset_img = (0.1, 0.1, 0.1)
    vec_points = [
        start_img_grid .- offset_img,
        finish_img .- offset_img,
        final_mid_cell_point .- offset_img,
    ]

    vec_intensities_hand = [intensity_array[1, 1, 1], intensity_array[end, end, end], hand_intensity]
    @test _eval_intensity(vec_points, f_img, offset=offset_img) ≈ f_img(vec_points, offset=offset_img) atol = 1e-1
    @test f_img(vec_points, offset=offset_img) ≈ vec_intensities_hand atol = 1e-1

end
