####################
# VTK Images Tests #
####################
using Test: @testset, @test

using Apolo.Images
using Apolo.Images.VTKImages
using Apolo.Images: _total_num_pixels, _extrapolate, _interpolate, _index_is_inbounds

using Apolo: vtk_structured_write, load_vtk_img

const INTERVAL_START = LinRange(-10.0, 10.0, 20)
const INTERVAL_POS = LinRange(1.0, 10.0, 20)
const INTERVAL_OFFSET = LinRange(-10.0, -1.0, 20)
const TOLERANCE = 1e-3
const DCM_TO_LOAD = "./test/DICOMImages/"


@testset "Images.VTKImage 2D" begin

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
        intensity_array, vtk_img_grid, spacing_img, start_img, path_img
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

@testset "Images.VTKImage 3D" begin

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
        intensity_array, vtk_img_grid, spacing_img, start_img, path_img
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
#     @test mtime == elapsed_time(img_data)
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
