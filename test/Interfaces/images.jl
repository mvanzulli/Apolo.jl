################
# Images tests #
################

# Using internal packages to test 
using Apolo.Images
using Apolo.ForwardProblem: getdim


# Using external packages to test 
using Test: @test, @testset
using HistogramThresholding: build_histogram
using Statistics:  var, mean
using DICOM: @tag_str


# Constant variables and tpyes
# this file is executed from /test/
const DCM_TO_LOAD = "./DICOMImages/"


# this tests is also validated with 3D slicer results
@testset "Reading DICOM" begin
    # load dicom images
    dcm = load_dicom(DCM_TO_LOAD)

    # get image dimensions
    dim = getdim(dcm);
    dims = (256, 256, 400)
    named_dims = (sagital = dims[1], coronal = dims[2], axial = dims[3])
    @test dim == named_dims

    # get image intensity array
    intensity = getintensity(dcm)
    @test size(intensity) == dims

    # get the image spacing
    spacing = (0.78125, 0.78125, 0.500003) 
    named_spacing = (sagital = spacing[1], coronal = spacing[2], axial = spacing[3])
    spacing = getspacing(dcm)
    @test spacing == named_spacing

        
    # test image orientation
    image_orientation = getorientation(dcm)
    image_orientation_vector = if image_orientation == :axial 
         image_orient = [1, 0, 0, 0, 1, 0]
    end

    # get all the image attributes and hyper parameters
    hyp = gethyp(dcm)
    patient_name = "COMPRESION MAMOGRAFICA PACIENTE 1 "
    @test hyp[tag"PatientName"] == patient_name == get_patient_name(dcm)
    @test hyp[tag"PixelSpacing"] == [named_spacing.sagital, named_spacing.coronal] 
    @test hyp[tag"SliceThickness"] == named_spacing.axial
    # test first slide with the last slice of the image parameters
    hyp[tag"PixelData"] == intensity[:,:,1]
    hyp_image_orientation = hyp[tag"ImageOrientationPatient"]
    @test ≈(hyp_image_orientation,  image_orientation_vector, atol = 1e-2)

end

@testset "Processing DICOM" begin

    # load dicom images
    dcm = load_dicom(DCM_TO_LOAD)

    intensity = getintensity(dcm)
    # build histogram using HistogramThresholding pkg
    num_segments = 3
    segments_pkg, num_pix_pkg = build_histogram(vec(intensity), num_segments)   
    # trasnform it into vectors
    segments_pkg = segments_pkg |> collect
    # remove the first elements 
    # compute segments 
    segments = compute_histogram(dcm, segments_pkg[2:end])
    # get segments indexes 
    seg_indexs = getindexes.(segments) 
    # get num pixels inside each segment
    num_pixes = getnumpix.(segments)
    # test results (offset array)
    @test isapprox(num_pix_pkg[1:end], num_pixes, atol=2)
    # get segments intensity vector 
    intensity_segments = getintensity.(segments)
    # compute means for each segments  
    mean_segments = mean.(segments)
    # compute std for each segments  
    std_segments = std.(segments)

end