################
# Images tests #
################
using Test: @testset, @test
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
