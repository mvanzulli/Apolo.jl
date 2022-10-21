################
# Images tests #
################

# External module functions 
using Apolo:
    create_rectangular_grid,
    dimension,
    getgrid,
    ⊂
# Internal exported module functions 
using Apolo.Images
# Internal non-exported module functions 
using Apolo.Images:
    _eval_intensity,
    _eval_intensity_inside_grid,
    _is_inside_img_grid,
    _index_is_inbounds,
    _node_cartesian_index,
    _build_ferrite_img_intensity,
    _point_is_inbounds,
    _which_border,
    _eval_intensity_outside_grid


# Using external packages to test 
using Test: @test, @testset
using AutoHashEquals
using HistogramThresholding: build_histogram
using Statistics: var, mean
using DICOM: @tag_str


# Define tolerances
const TOLERANCE = 1e-3

@testset "Generic images unitary tests" begin

    # TODO: add random values
    # 2D image 
    # image dimensions 
    spacing_img = (0.5, 0.25)
    num_pixels = (4, 3)
    start_img = (1.0, 1.0)
    image_dimension = length(spacing_img)
    finish_img = start_img .+ num_pixels .* spacing_img
    length_img = finish_img .- start_img
    # image intensity
    intensity_vec = collect(1.0:Float64(prod(num_pixels)))
    intensity_mat = reshape(intensity_vec, num_pixels)
    # create the image grid
    grid_img = create_img_grid(start_img, spacing_img, length_img, num_pixels)

    # generate image 
    generic_img = GenericImage(intensity_mat, num_pixels, spacing_img, start_img)

    # test extractors
    @test intensity(generic_img) == intensity_mat
    @test intensity_type(generic_img) == eltype(intensity_mat)
    @test dimension(generic_img) == image_dimension
    @test numpix(generic_img) == num_pixels
    @test total_numpix(generic_img) == prod(num_pixels)
    @test spacing(generic_img) == spacing_img
    @test start(generic_img) == start_img
    @test finish(generic_img) == finish_img
    @test length(generic_img) == length_img

    # test Checks
    # indexes 
    # is inside function
    @test !(_index_is_inbounds(-1, generic_img))
    @test !(_index_is_inbounds(total_numpix(generic_img) + 1, generic_img))
    # normal case 
    @test _node_cartesian_index(10, generic_img) == CartesianIndex(2, 3)
    # edge cases
    @test _node_cartesian_index(2, generic_img) == CartesianIndex(2, 1)
    @test _node_cartesian_index(12, generic_img) == CartesianIndex(4, 3)
    # points
    @test _point_is_inbounds(start_img, generic_img)
    @test _point_is_inbounds(finish_img, generic_img)
    @test _point_is_inbounds((start_img .+ finish_img) ./ 2, generic_img)
    @test !(_point_is_inbounds((start_img .+ (1.0, 1.0)), generic_img))

    # test coordinates generator
    coords = coordinates(generic_img)
    @test spacing_img[1] == coords[1][2] - coords[1][1]
    @test spacing_img[2] == coords[2][2] - coords[2][1]
    extremaᵢ, extremaⱼ = extrema.(coords)
    @test extremaᵢ == (start_img[1] + spacing_img[1] / 2, finish_img[1] - spacing_img[1] / 2)
    @test extremaⱼ == (start_img[2] + spacing_img[2] / 2, finish_img[2] - spacing_img[2] / 2)

    # test mesh generator
    grid_img_extracted = getgrid(generic_img)
    @test grid_img.cells == grid_img_extracted.cells
    @test grid_img.nodes == grid_img_extracted.nodes
    @test grid_img.facesets == grid_img_extracted.facesets

    # check grid image methods
    inside_grid_img_point = (1.5, 1.5)

    # build ferrite vector to evaluate via PointEvalHandler
    _build_ferrite_img_intensity(generic_img, grid_img)

    # Evaluate intensity via an image functor 
    # p is inside the grid (interpolate)
    @test inside_grid_img_point ⊂ grid_img && _is_inside_img_grid(generic_img, inside_grid_img_point)
    img_grid_origin = start_img .+ spacing_img ./ 2
    @test start_grid(generic_img) == img_grid_origin
    img_grid_finish = finish_img .- spacing_img ./ 2
    @test finish_grid(generic_img) == img_grid_finish

    # test origin intensity
    @test (coords[1][1], coords[2][1]) == img_grid_origin
    @test _eval_intensity(generic_img, img_grid_origin)[1] ==
          generic_img(img_grid_origin...)[1] ==
          _eval_intensity_inside_grid(generic_img, img_grid_origin)[1] == intensity_vec[1]

    # test second point (start_grid + (Δₓ/2, 0)) intensity
    consecutive_x_in = img_grid_origin .+ (spacing_img[1], 0) ./ 2
    @test _eval_intensity(generic_img, consecutive_x_in) ==
          generic_img(consecutive_x_in...) ==
          _eval_intensity_inside_grid(generic_img, consecutive_x_in) ==
          [(intensity_vec[1] + intensity_vec[2]) / 2]

    # test second point (start_grid + (0, Δⱼ/2)) intensity
    consecutive_y_in = img_grid_origin .+ (0, spacing_img[2]) ./ 2
    @test _eval_intensity(generic_img, consecutive_y_in) ==
          generic_img(consecutive_y_in...) ==
          _eval_intensity_inside_grid(generic_img, consecutive_y_in) == [intensity_vec[1] + intensity_vec[5]] / 2

    # at the middle of the end cell (finish - (Δᵢ, Δⱼ)/2)
    mid_final_cell = img_grid_finish .- spacing_img ./ 2
    intensity_mid_end_cell = (intensity_vec[7] + intensity_vec[8] + intensity_vec[11] + intensity_vec[12]) / 4
    @test _eval_intensity(generic_img, mid_final_cell) ==
          generic_img(mid_final_cell...) ==
          _eval_intensity_inside_grid(generic_img, mid_final_cell) ==
          [intensity_mid_end_cell]

    # TODO: p is not in the grid (extrapolate)
    # Test corners detection 
    p_left_bottom = start(generic_img)
    p_right_top = finish(generic_img)
    p_left_top = (start(generic_img)[1], finish(generic_img)[2])
    p_right_bottom = (finish(generic_img)[1], start(generic_img)[2])
    
    @test _which_border(generic_img, p_left_bottom) == :left_bottom
    @test _which_border(generic_img, p_right_top) == :right_top
    @test _which_border(generic_img, p_left_top) == :left_top
    @test _which_border(generic_img, p_right_bottom) == :right_bottom
    div = rand(2:10)
    xₚ = start(generic_img)[1] + ((start(generic_img)[1] + finish(generic_img)[1]) ./ div)
    yₚ = start(generic_img)[2] + ((start(generic_img)[2] + finish(generic_img)[2]) ./ div)
    
    # Test borders detection 
    p_bottom = (xₚ, start(generic_img)[2])
    p_top = (xₚ, finish(generic_img)[2])
    p_left = (start(generic_img)[1], yₚ)
    p_right = (finish(generic_img)[1], yₚ)
    @test _which_border(generic_img, p_bottom) == :bottom
    @test _which_border(generic_img, p_top) == :top
    @test _which_border(generic_img, p_left) == :left
    @test _which_border(generic_img, p_right) == :right

    @test _eval_intensity_outside_grid(generic_img, p_left_bottom) == intensity(generic_img)[1, 1]
    @test _eval_intensity_outside_grid(generic_img, p_right_top) == intensity(generic_img)[end, end]
    @test _eval_intensity_outside_grid(generic_img, p_left_top) == intensity(generic_img)[1, end]
    @test _eval_intensity_outside_grid(generic_img, p_right_bottom) == intensity(generic_img)[end, 1]

    # Test bottom outside interpolation
    # test point (point_bottom_x + n(Δᵢ, 0)) intensity
    batch_x = rand(1:num_pixels[1] -2)

    # point that matches with a node
    point_bottom_x = (p_left_bottom .+ batch_x .* (spacing_img[1], 0))
    @test position = _which_border(generic_img, point_bottom_x) == :bottom
    @test _eval_intensity_outside_grid(generic_img, point_bottom_x) ==
    generic_img(point_bottom_x...)[1] ==
    intensity(generic_img)[batch_x + 1, 1]

    # point that doesn't matches with a node
    fraction = 2/3
    point_bottom_x_inter = (point_bottom_x .+ batch_x .* (spacing_img[1], 0) .* fraction)
    @test position = _which_border(generic_img, point_bottom_x_inter) == :bottom
    @test _eval_intensity_outside_grid(generic_img, point_bottom_x_inter) ==
    generic_img(point_bottom_x_inter...)[1]
    hand_interpolation = intensity(generic_img)[batch_x + 2, 1] * fraction + intensity(generic_img)[batch_x + 1, 1] *(1-fraction)
    @test isapprox(hand_interpolation, generic_img(point_bottom_x_inter...)[1], atol = TOLERANCE )
    
    # Test top outside interpolation
    # test point (p_top_x + n(Δᵢ, end)) intensity
    batch_x = rand(1:num_pixels[1] -2)

    # point that matches with a node
    point_top_x = (p_top_x .+ batch_x .* (spacing_img[1], 0))
    @test position = _which_border(generic_img, point_top_x) == :top
    @test _eval_intensity_outside_grid(generic_img, point_top_x) ==
    generic_img(point_top_x...)[1] ==
    intensity(generic_img)[batch_x + 1, end]

    # point that doesn't matches with a node
    point_top_x_inter = (point_top_x .+ batch_x .* (spacing_img[1], 0) .* fraction)
    @test position = _which_border(generic_img, point_top_x_inter) == :top
    @test _eval_intensity_outside_grid(generic_img, point_top_x_inter) ==
    generic_img(point_top_x_inter...)[1]
    hand_interpolation = intensity(generic_img)[batch_x + 2, end] * fraction + intensity(generic_img)[batch_x + 1, end] *(1-fraction)
    @test isapprox(hand_interpolation, generic_img(point_top_x_inter...)[1], atol = TOLERANCE )


    # Test left outside interpolation
    # test point (p_left_bottom + n(0, Δⱼ)) intensity
    batch_y = rand(1:num_pixels[2] -2)

    # point that matches with a node
    point_left_y = (p_left_bottom .+ batch_y .* (0, spacing_img[2]))
    position = _which_border(generic_img, point_left_y) == :left
    @test _eval_intensity_outside_grid(generic_img, point_left_y) ==
    generic_img(point_left_y...)[1] ==
    intensity(generic_img)[1, batch_y + 1]
    
    # point that doesn't matches with a node
    point_left_y_inter = (point_left_y .+ batch_y .* (0, spacing_img[2]) .* fraction)
    @test position = _which_border(generic_img, point_left_y_inter) == :left
    @test _eval_intensity_outside_grid(generic_img, point_left_y_inter) ==
    generic_img(point_left_y_inter...)[1]
    hand_interpolation = intensity(generic_img)[1, batch_y + 2] * fraction + intensity(generic_img)[1, batch_y + 1] *(1-fraction)
    @test isapprox(hand_interpolation, generic_img(point_left_y_inter...)[1], atol = TOLERANCE )
    
    # Test right outside interpolation
    # test point (p_right_bottom + n(0, Δⱼ)) intensity
    batch_y = rand(1:num_pixels[2] -1)
    # point that matches with a node
    point_right_y = (p_right_bottom .+ batch_y .* (0, spacing_img[2]))
    position = _which_border(generic_img, point_right_y) == :right
    @test _eval_intensity_outside_grid(generic_img, point_right_y) ==
    generic_img(point_right_y...)[1] ==
    intensity(generic_img)[end, batch_y + 1]
    
    # point that doesn't matches with a node
    point_right_y_inter = (point_right_y .+ batch_y .* (0, spacing_img[2]).* fraction)
    position = _which_border(generic_img, point_right_y_inter) == :right
    hand_interpolation = intensity(generic_img)[end, batch_y + 2] * fraction + intensity(generic_img)[end, batch_y + 1] *(1-fraction)
    @test _eval_intensity_outside_grid(generic_img, point_right_y_inter) ==
    generic_img(point_right_y_inter...)[1]
    @test isapprox(hand_interpolation, generic_img(point_right_y_inter...)[1], atol = TOLERANCE )

end


# Constant variables and tpyes
# this file is executed from /test/
const DCM_TO_LOAD = "./DICOMImages/"


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


