################
# Images tests #
################
using Apolo, Test
using Apolo.Geometry: _convert_to_ferrite_nomenclature
using Apolo.Images
using Apolo.Images: _total_num_pixels, _index_is_inbounds, _eval_intensity,
    _interpolate, _extrapolate

# using AutoHashEquals
# using HistogramThresholding: build_histogram
# using Statistics: var, mean
# using DICOM: @tag_str

const INTERVAL_START = LinRange(-10.0, 10.0, 20)
const INTERVAL_POS = LinRange(1.0, 10.0, 20)
const INTERVAL_OFFSET = LinRange(-10.0, -1.0, 20)
const TOLERANCE = 1e-3
const DCM_TO_LOAD = "./test/DICOMImages/"


@testset "Ferrite 2D image unitary tests" begin

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

    # generate defualt image
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

@testset "Ferrite 3D image unitary tests" begin

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

    # generate image and internaly cash a grid
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

@testset "Analytic 3D image unitary tests" begin

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

@testset "VTK 2D image unitary tests" begin

    # Define VTK image properties
    intensity_function(x, y) = 2x + 2y

    start_img = (0.0, 0.0)
    spacing_img = (0.25, 0.25)
    num_pixels_img = (4, 4)
    start_img_grid = start_img .+ spacing_img ./ 2
    image_dimension = length(spacing_img)
    finish_img = start_img .+ num_pixels_img .* spacing_img
    finish_img_grid = finish_img .- spacing_img ./ 2
    length_img = finish_img .- start_img


    xc = LinRange(start_img_grid[1], finish_img_grid[1], num_pixels_img[1])
    yc = LinRange(start_img_grid[2], finish_img_grid[2], num_pixels_img[2])
    coords = [xc, yc]

    intensity_array = [intensity_function(x, y) for x in xc for y in yc]

    # Write a VTK image
    path_img = tempname()
    @info "Writing VTK 2D image in $path_img"
    vtk_structured_write(coords, intensity_function, :intensity, path_img, "")
    intensity_array = reshape(intensity_array, num_pixels_img)
    vtk_structured_write(coords, intensity_array, :intensity, path_img, "")

    # Test VTK image features
    vtk_img = VTKImage(
        intensity_array, spacing_img, start_img, path_img, ferrite_grid=true
    )
    @test path(vtk_img) == path_img

    # Build a VTK image with a given grid
    vtk_img_grid = grid(vtk_img)
    vtk_img = VTKImage(
        intensity_array, vtk_img_grid,  spacing_img, start_img, path_img
        )
    # Read a VTK image
    vtk_img_red = load_vtk_img(path_img)

    # Test written and red images are the same
    @test intensity(vtk_img_red) == intensity(vtk_img)
    @test num_pixels(vtk_img_red) == num_pixels(vtk_img)
    @test finish(vtk_img_red) == finish(vtk_img)
    @test finish_grid(vtk_img_red) == finish_grid(vtk_img)
    @test start(vtk_img_red) == start(vtk_img)
    @test start_grid(vtk_img_red) == start_grid(vtk_img)
    @test spacing(vtk_img_red) == spacing(vtk_img)
    @test coordinates(vtk_img_red) == coordinates(vtk_img)

    # test getter functions
    @test start(vtk_img) == start_img
    @test start_grid(vtk_img) == start(grid(vtk_img)) == start_img_grid
    @test finish(vtk_img) == finish_img
    @test finish_grid(vtk_img) == finish(grid(vtk_img)) == finish_img_grid
    @test collect(length(vtk_img)) ≈ collect(length_img) atol = TOLERANCE
    coords = coordinates(vtk_img)
    @test coords[1] == LinRange.(start_img_grid, finish_img_grid, num_pixels_img)[1]
    @test coords[2] == LinRange.(start_img_grid, finish_img_grid, num_pixels_img)[2]
    @test spacing_img[1] ≈ (coords[1][2] - coords[1][1]) atol = TOLERANCE
    @test spacing_img[2] ≈ (coords[2][2] - coords[2][1]) atol = TOLERANCE

    @test extrema(vtk_img) == [
        (start_img[1], finish_img[1]),
        (start_img[2], finish_img[2]),
    ]

    @test intensity(vtk_img) == intensity_array
    @test intensity_type(vtk_img) == eltype(value(vtk_img.intensity))
    @test dimension(vtk_img) == image_dimension
    @test num_pixels(vtk_img) == num_pixels_img
    @test _total_num_pixels(vtk_img) == prod(num_pixels_img)
    @test spacing(vtk_img) == spacing_img

    # test checker functions
    # ∈ and ∉
    @test start_img ∈ vtk_img
    @test finish_img ∈ vtk_img
    @test (start_img .+ (finish_img .- start_img) ./ 2) ∈ vtk_img
    @test (start_img .+ 900 .* finish_img) ∉ vtk_img

    # final point the image grid
    @test _extrapolate([finish_img .+ (1.0, 1.0)], vtk_img) ≈ [intensity_array[end, end]] atol = TOLERANCE
    @test [intensity_array[end, end]] ≈ [vtk_img(finish_img...)] atol = TOLERANCE

    # test point at the final cell middle point
    final_mid_cell_point = finish_img_grid .- spacing_img ./ 2
    hand_intensity = sum(intensity_array[end-1:end, end-1:end]) / 4
    @test _interpolate([final_mid_cell_point], vtk_img) ≈ [hand_intensity] atol = TOLERANCE
    @test [hand_intensity] ≈ [vtk_img(final_mid_cell_point...)] atol = TOLERANCE

    # test with a random point
    rand_point = start_img .+ length_img ./ rand(1:10)
    @test vtk_img(rand_point...) ≈ intensity_function(rand_point...) atol = TOLERANCE skip = true

end

@testset "VTK 3D image unitary tests" begin

    # Define VTK image properties
    intensity_function(x, y, z) = 2x - 3y + 4z

    start_img = (rand(INTERVAL_START), rand(INTERVAL_START), rand(INTERVAL_START))
    spacing_img = (rand(INTERVAL_POS), rand(INTERVAL_POS), rand(INTERVAL_POS)) ./ 100
    num_pixels_img = (15, 15, 15)
    start_img_grid = start_img .+ spacing_img ./ 2
    image_dimension = length(spacing_img)
    finish_img = start_img .+ num_pixels_img .* spacing_img
    finish_img_grid = finish_img .- spacing_img ./ 2
    length_img = finish_img .- start_img


    xc = LinRange(start_img_grid[1], finish_img_grid[1], num_pixels_img[1])
    yc = LinRange(start_img_grid[2], finish_img_grid[2], num_pixels_img[2])
    zc = LinRange(start_img_grid[3], finish_img_grid[3], num_pixels_img[3])
    coords = [xc, yc, zc]

    intensity_array = [intensity_function(x, y, z) for x in xc for y in yc for z in zc]
    intensity_array = reshape(intensity_array, num_pixels_img)

    # Write a VTK image (structured grids only)
    path_img = tempname()
    @info "Writing VTK 3D image in $path_img"
    vtk_structured_write(coords, intensity_function, :intensity, path_img, "")
    vtk_structured_write(coords, intensity_array, :intensity, path_img, "")

    # Test VTK image features
    vtk_img = VTKImage(
        intensity_array, spacing_img, start_img, path_img, ferrite_grid=true
    )
    # Build a VTK image with a given grid
    vtk_img_grid = grid(vtk_img)
    vtk_img = VTKImage(
        intensity_array, vtk_img_grid,  spacing_img, start_img, path_img
        )

    # Read a VTK image
    vtk_img_red = load_vtk_img(path_img)

    # Test written and red images are the same
    @test intensity(vtk_img_red) ≈ intensity(vtk_img) atol = TOLERANCE
    @test num_pixels(vtk_img_red) == num_pixels(vtk_img)
    @test collect(finish(vtk_img_red)) ≈ collect(finish(vtk_img)) atol = TOLERANCE
    @test collect(finish_grid(vtk_img_red)) ≈ collect(finish_grid(vtk_img)) atol = TOLERANCE
    @test collect(start(vtk_img_red)) ≈ collect(start(vtk_img)) atol = TOLERANCE
    @test collect(start_grid(vtk_img_red)) ≈ collect(start_grid(vtk_img)) atol = TOLERANCE
    @test collect(spacing(vtk_img_red)) ≈ collect(spacing(vtk_img)) atol = TOLERANCE

    # test getter functions
    @test start(vtk_img) == start_img
    @test start_grid(vtk_img) == start(grid(vtk_img)) == start_img_grid
    @test finish(vtk_img) == finish_img
    @test finish_grid(vtk_img) == finish(grid(vtk_img)) == finish_img_grid
    @test collect(length(vtk_img)) ≈ collect(length_img) atol = TOLERANCE
    coords = coordinates(vtk_img)
    @test coords[1] == LinRange.(start_img_grid, finish_img_grid, num_pixels_img)[1]
    @test coords[2] == LinRange.(start_img_grid, finish_img_grid, num_pixels_img)[2]
    @test coords[3] == LinRange.(start_img_grid, finish_img_grid, num_pixels_img)[3]
    @test spacing_img[1] ≈ (coords[1][2] - coords[1][1]) atol = TOLERANCE
    @test spacing_img[2] ≈ (coords[2][2] - coords[2][1]) atol = TOLERANCE
    @test spacing_img[3] ≈ (coords[3][2] - coords[3][1]) atol = TOLERANCE

    @test extrema(vtk_img) == [
        (start_img[1], finish_img[1]),
        (start_img[2], finish_img[2]),
        (start_img[3], finish_img[3])
    ]

    @test intensity(vtk_img) == intensity_array
    @test intensity_type(vtk_img) == eltype(value(vtk_img.intensity))
    @test dimension(vtk_img) == image_dimension
    @test num_pixels(vtk_img) == num_pixels_img
    @test _total_num_pixels(vtk_img) == prod(num_pixels_img)
    @test spacing(vtk_img) == spacing_img

    # test checker functions
    # ∈ and ∉
    @test start_img ∈ vtk_img
    @test finish_img ∈ vtk_img
    @test (start_img .+ (finish_img .- start_img) ./ 2) ∈ vtk_img
    @test (start_img .+ 900 .* finish_img) ∉ vtk_img

    # indexes
    # is inside function
    @test !(_index_is_inbounds(-1, vtk_img))
    @test !(_index_is_inbounds(_total_num_pixels(vtk_img) + 1, vtk_img))

    @test _interpolate([start_img_grid], vtk_img) ≈ [intensity_array[1, 1, 1]] atol = TOLERANCE
    @test [vtk_img(start_img_grid...)] ≈ [intensity_array[1, 1, 1]] atol = TOLERANCE

    # final point the image grid
    @test _extrapolate([finish_img .+ (1.0, 1.0, 1.0)], vtk_img) ≈ [intensity_array[end, end, end]] atol = TOLERANCE
    @test [intensity_array[end, end, end]] ≈ [vtk_img(finish_img...)] atol = TOLERANCE

    # test point at the final cell middle point
    final_mid_cell_point = finish_img_grid .- spacing_img ./ 2
    hand_intensity = sum(intensity_array[end-1:end, end-1:end, end-1:end]) / 8
    @test _interpolate([final_mid_cell_point], vtk_img) ≈ [hand_intensity] atol = TOLERANCE
    @test [hand_intensity] ≈ [vtk_img(final_mid_cell_point...)] atol = TOLERANCE

    # test with a random point
    rand_point = start_img .+ length_img ./ rand(1:100)
    itp_rand_int = vtk_img(rand_point...)
    exact_rand_int = intensity_function(rand_point...)
    @test exact_rand_int ≈ itp_rand_int atol = abs(itp_rand_int * 1e-1)

end

# @testset "VTK 3D sequence" begin

#     # Define VTK image properties
#     intensity_function(x, y, z, t) = 2*sin(t) * cos(9x - 1y + 7z)

#     start_img = (rand(INTERVAL_START), rand(INTERVAL_START), rand(INTERVAL_START))
#     spacing_img = (rand(INTERVAL_POS), rand(INTERVAL_POS), rand(INTERVAL_POS)) ./ 100
#     num_pixels_img = (8, 5, 11)
#     start_img_grid = start_img .+ spacing_img ./ 2
#     image_dimension = length(spacing_img)
#     finish_img = start_img .+ num_pixels_img .* spacing_img
#     finish_img_grid = finish_img .- spacing_img ./ 2
#     length_img = finish_img .- start_img
#     mtime = LinRange(rand(INTERVAL_POS), 3*rand(INTERVAL_POS) , 8)

#     xc = LinRange(start_img_grid[1], finish_img_grid[1], num_pixels_img[1])
#     yc = LinRange(start_img_grid[2], finish_img_grid[2], num_pixels_img[2])
#     zc = LinRange(start_img_grid[3], finish_img_grid[3], num_pixels_img[3])
#     coords = [xc,yc,zc]
#     vars = [coords..., mtime]

#     intensity_array = [intensity_function(var...) for var in Iterators.product(vars...)]

#     # Write the VTK sequence of images
#     tdir = joinpath(tempdir(),"temp")
#     mkdir(tdir)
#     tname = splitpath(tempname())[end]
#     tempfile = joinpath(tdir, tname)
#     root_name_sequence = "temp"

#     @info "Writing VTK 3D sequence in $tempdir"
#     vtk_structured_write_sequence(coords, intensity_array, :intensity, tname, tdir)
#     # vtk_structured_write_sequence(vars, intensity_function, :intensity, tname, tdir)

#     # Read the VTK sequence in a vector of images
#     imgs = load_vtk_sequence_imgs(tdir)
#     rm(tdir, recursive = true, force = true)

#     # Check that every image have the same grid on memory
#     rand_t_idx = rand(2:length(mtime))
#     rand_img = imgs[rand_t_idx]
#     @test grid(imgs[1]) === grid(rand_img)
#     rand_point_1 = start_img_grid .+ (length_img ./ rand(2:10))
#     int_rand_1 = intensity_function(rand_point_1...,mtime[rand_t_idx])
#     rand_point_2 = start_img_grid .+ (length_img ./ rand(2:10))
#     int_rand_2 = intensity_function(rand_point_2...,mtime[rand_t_idx])

#     ints_loaded = rand_img([rand_point_1, rand_point_2])
#     @test [int_rand_1, int_rand_2] ≈ ints_loaded atol = .5

#     # Define an Image data and test
#     start_roi = start_img_grid
#     finish_roi = (start_img_grid .+ length_img ./ 4)
#     roi_func(x) = all(start_roi .≤ (x[1], x[2], x[3]) .≤ finish_roi)
#     img_data = ImageData(imgs, roi_func, mtime )
#     img_ref = reference_img(img_data)
#     imgs_def = deformed_imgs(img_data)
#     rand_img_def = imgs_def[rand_t_idx-1]
#     @test roi_func == roi(img_data)
#     @test mtime == time_measured(img_data)
#     @test grid(img_data) === grid(img_ref) === grid(rand_img_def)

#     # Get roi coorindates
#     nodes_roi = roi_nodes(img_data)
#     @test collect(nodes_roi) == collect(grid(grid(img_data)).nodesets["roi"])
#     roi_vec_coords = roi_nodes_coords(img_data)
#     @test all([roi_func(node_roi_coord) for node_roi_coord in roi_vec_coords])

#     # Get roi intensity and check values reference image
#     int_ref_roi_to_test = img_ref(roi_vec_coords)
#     int_ref_roi_bench = [intensity_function(p...,mtime[1]) for p in roi_vec_coords]
#     @test int_ref_roi_to_test ≈ int_ref_roi_bench atol = TOLERANCE

#     # Get roi intensity and check values random defomed image
#     int_def_roi_to_test = rand_img_def(roi_vec_coords)
#     int_def_roi_bench = [intensity_function(p...,mtime[rand_t_idx]) for p in roi_vec_coords]
#     @test int_ref_roi_to_test ≈ int_ref_roi_bench atol = TOLERANCE
# end

@testset "DICOM 3D image unitary tests" begin

    # Define VTK image properties
    intensity_function(x, y, z) = 2x - 3y + 4z

    start_img = (rand(INTERVAL_START), rand(INTERVAL_START), rand(INTERVAL_START))
    spacing_img = (rand(INTERVAL_POS), rand(INTERVAL_POS), rand(INTERVAL_POS)) ./ 50
    num_pixels_img = (axial=4, sagital=3, radial=2)
    start_img_grid = start_img .+ spacing_img ./ 2
    image_dimension = length(spacing_img)
    finish_img = Tuple(start_img .+ collect(num_pixels_img) .* spacing_img)
    finish_img_grid = finish_img .- spacing_img ./ 2
    length_img = finish_img .- start_img


    xc = LinRange(start_img_grid[1], finish_img_grid[1], collect(num_pixels_img)[1])
    yc = LinRange(start_img_grid[2], finish_img_grid[2], collect(num_pixels_img)[2])
    zc = LinRange(start_img_grid[3], finish_img_grid[3], collect(num_pixels_img)[3])
    coords = [xc, yc, zc]

    intensity_array = [intensity_function(x, y, z) for x in xc for y in yc for z in zc]
    intensity_array = reshape(intensity_array, Tuple(collect(num_pixels_img)))

    # Test constructor
    med_img = MedicalImage(intensity_array, spacing_img, start_img)

    # Special medical getter functions
    @test hyper_parameters(med_img) == nothing
    @test path(med_img) == ""
    @test orientation(med_img) == [:sagital, :coronal, :axial]

    # test getter functions
    @test start(med_img) == start_img
    @test start_grid(med_img) == start(grid(med_img)) == start_img_grid
    @test finish(med_img) == finish_img
    @test finish_grid(med_img) == finish(grid(med_img)) == finish_img_grid
    @test collect(length(med_img)) ≈ collect(length_img) atol = TOLERANCE
    coords = coordinates(med_img)
    @test spacing_img[1] ≈ (coords[1][2] - coords[1][1]) atol = TOLERANCE
    @test spacing_img[2] ≈ (coords[2][2] - coords[2][1]) atol = TOLERANCE
    @test spacing_img[3] ≈ (coords[3][2] - coords[3][1]) atol = TOLERANCE

    @test extrema(med_img) == [
        (start_img[1], finish_img[1]),
        (start_img[2], finish_img[2]),
        (start_img[3], finish_img[3])
    ]
    @test intensity(med_img) == intensity_array
    @test intensity_type(med_img) == intensity_type(med_img.intensity)
    @test dimension(med_img) == image_dimension
    @test num_pixels(med_img) == Tuple(collect(num_pixels_img))
    @test _total_num_pixels(med_img) == prod(num_pixels_img)
    @test spacing(med_img) == spacing_img

    # test checker functions
    # ∈ and ∉
    @test start_img ∈ med_img
    @test finish_img ∈ med_img
    @test (start_img .+ (finish_img .- start_img) ./ 2) ∈ med_img
    @test (start_img .+ 9 .* finish_img) ∉ med_img

    # FIXME: IMPROVE PERFORMANCE
    # # Load a .dcm folder
    # time = @elapsed begin
    #     med_img = load_medical_img(DCM_TO_LOAD)
    # end

end


#=

# Constant variables and tpyes
# this file is executed from /test/


# this tests is also validated with 3D slicer results
@testset "Reading DICOM" begin
    # load dicom images
    med_img = load_dicom(DCM_TO_LOAD)

    # get image dimensions
    dim = dimension(med_img)
    num_pix = (256, 256, 400)
    named_num_pix = (sagital=num_pix[1], coronal=num_pix[2], axial=num_pix[3])
    @test dim == named_num_pix

    # get image intensity array
    intensity_med = intensity(med_img)
    @test size(intensity_med) == num_pix

    # get the image spacing
    spacing_test = (0.78125, 0.78125, 0.500003)
    named_spacing = (
        sagital=spacing_test[1],
        coronal=spacing_test[2],
        axial=spacing_test[3]
    )
    spacing_med = spacing(med_img)
    @test spacing_med == named_spacing


    # test image orientation
    image_orientation = orientation(med_img)
    image_orientation_vector = if image_orientation == :axial
        image_orient = [1, 0, 0, 0, 1, 0]
    end

    # get all the image attributes and hyper parameters
    hyp = hyper_parameters(med_img)
    patient_name_test = "COMPRESION MAMOGRAFICA PACIENTE 1 "
    @test hyp[tag"PatientName"] == patient_name_test == patient_name(med_img)
    @test hyp[tag"PixelSpacing"] == [named_spacing.sagital, named_spacing.coronal]
    @test hyp[tag"SliceThickness"] == named_spacing.axial
    # test first slide with the last slice of the image parameters
    hyp[tag"PixelData"] == intensity_med[:, :, 1]
    hyp_image_orientation = hyp[tag"ImageOrientationPatient"]
    @test ≈(hyp_image_orientation, image_orientation_vector, atol=1e-2)

end

@testset "Processing DICOM" begin

    # load dicom images
    med_img = load_dicom(DCM_TO_LOAD)

    intensity_med = intensity(med_img)
    # build histogram using HistogramThresholding pkg
    num_segments = 3
    segments_pkg, num_pix_pkg = build_histogram(vec(intensity_med), num_segments)
    # transform it into vectors
    segments_pkg = segments_pkg |> collect
    # remove the first elements
    # compute segments
    segments = compute_histogram(med_img, segments_pkg[2:end])
    # get segments indexes
    seg_indexs = indexes.(segments)
    # get num pixels inside each segment
    num_pixes = numpix.(segments)
    # test results (offset array)
    @test isapprox(num_pix_pkg[1:end], num_pixes, atol=2)
    # get segments intensity vector
    intensity_segments = intensity.(segments)
    # compute means for each segments
    mean_segments = mean.(segments)
    # compute std for each segments
    std_segments = std.(segments)

end


=#
