################
# Images tests #
################
using Apolo.Images
using Apolo.Geometry

using Apolo.Geometry: _convert_to_ferrite_nomenclature
using Apolo.Images:_total_num_pixels, _index_is_inbounds, _eval_intensity, _interpolate, _extrapolate


# Using external packages to test
using Test: @test, @testset



# using AutoHashEquals
# using HistogramThresholding: build_histogram
# using Statistics: var, mean
# using DICOM: @tag_str


# Define tolerances

const INTERVAL_START = LinRange(-10.0, 10.0, 20)
const TOLERANCE = 1e-3

@testset "Ferrite images unitary tests" begin

    ################
    # 2D image test
    ##################
    # image dimensions
    start_img = (rand(INTERVAL_START), rand(INTERVAL_START))
    spacing_img = (0.5, 0.25)
    num_pixels_img = (4, 3)
    start_img_grid = start_img .+ spacing_img./2
    image_dimension = length(spacing_img)
    finish_img = start_img .+ num_pixels_img .* spacing_img
    finish_img_grid = finish_img .- spacing_img./2
    @test finish(start_img, num_pixels_img, spacing_img) == finish_img
    length_img = finish_img .- start_img
    @test length(start_img, finish_img) == length_img
    # image intensity
    intensity_vec = collect(1.0:Float64(prod(num_pixels_img)))
    intensity_array = reshape(intensity_vec, num_pixels_img)
    # create the image ferrite grid
    fgrid_img = create_ferrite_img_fgrid(start_img, spacing_img, length_img, num_pixels_img)
    intensity_ferrite_nomenclature = _convert_to_ferrite_nomenclature(intensity_array, fgrid_img)

    # generate image
    f_img = FerriteImage(intensity_array, num_pixels_img, start_img, spacing_img)
    @test grid(f_img).grid.cells == fgrid_img.grid.cells

    # test getter functions
    @test start(f_img) == start_img
    @test start_grid(f_img) == start(grid(f_img)) == start_img_grid
    @test finish(f_img) == finish_img
    @test finish_grid(f_img) == finish(grid(f_img)) == finish_img_grid
    @test length(f_img) == length_img
    coords = coordinates(f_img)
    @test coords == LinRange.(start_img_grid, finish_img_grid, num_pixels_img)
    @test isapprox(spacing_img[1], coords[1][2] - coords[1][1], atol = TOLERANCE )
    @test isapprox(spacing_img[2], coords[2][2] - coords[2][1], atol = TOLERANCE )

    @test extrema(f_img) == [(start_img[1], finish_img[1]), (start_img[2], finish_img[2])]
    @test intensity(f_img) == intensity_array
    @test intensity_type(f_img) == intensity_type(f_img.fintensity)
    @test dimension(f_img) == image_dimension
    @test num_pixels(f_img) == num_pixels_img
    @test _total_num_pixels(f_img) == prod(num_pixels_img)
    @test spacing(f_img) == spacing_img

    # test checker functions
    # ⊂ and ⊄
    @test start_img ⊂ f_img
    @test finish_img ⊂ f_img
    @test (start_img .+ (finish_img .- start_img) ./ 2) ⊂ f_img
    @test (start_img .+ 9 .* finish_img) ⊄ f_img

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
    @test isapprox(_interpolate([start_img_grid], f_img), [intensity_array[1,1]], atol = TOLERANCE )
    @test isapprox([f_img(start_img_grid...)], [intensity_array[1,1]], atol = TOLERANCE )

    # final poit the image grid
    @test _extrapolate([finish_img], f_img) == [intensity_array[end,end]] == [f_img(finish_img...)]

    # test point at the final cell middle point
    final_mid_cell_point = finish_img_grid .- spacing_img ./2
    hand_intensity = sum(intensity_array[end-1: end, end-1: end])/4
    @test _interpolate([final_mid_cell_point], f_img) == [hand_intensity] ==
    [f_img(final_mid_cell_point...)]

end


#=

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


=#
