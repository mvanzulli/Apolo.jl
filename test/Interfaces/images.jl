################
# Images tests #
################

# Using internal packages to test 
using Apolo.Images

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
    intensity_array = getintensity(dcm)
    @test size(intensity_array) == dims

    # get the image spacing
    spacing = (0.78125, 0.78125, 0.500003) 
    named_spacing = (sagital = spacing[1], coronal = spacing[2], axial = spacing[3])
    spacing = getspacing(dcm)
    @test spacing == named_spacing

    # get the image intensity array 
    int_array = getintensity(dcm)
        
    # test image orientation
    image_orientation = getorientation(dcm)
    image_orientation_vector = if image_orientation == :axial 
         image_orient = [1, 0, 0, 0, 1, 0]
    end

    # get all the image attributes and hyper parameters
    hyp = gethyp(dcm)
    patient_name = "COMPRESION MAMOGRAFICA PACIENTE 1 "
    @test hyp[tag"PatientName"] == patient_name
    @test hyp[tag"PixelSpacing"] == [named_spacing.sagital, named_spacing.coronal] 
    @test hyp[tag"SliceThickness"] == named_spacing.axial
    # test first slide with the last slice of the image parameters
    hyp[tag"PixelData"] == int_array[:,:,1]
    hyp_image_orientation = hyp[tag"ImageOrientationPatient"]
    @test ≈(hyp_image_orientation,  image_orientation_vector, atol = 1e-2)

end

@testset "Processing DICOM" begin

    # load dicom images
    dcm = load_dicom(DCM_TO_LOAD)

    intensity_array = getintensity(dcm)
    # build histogram using HistogramThresholding pkg
    num_segments = 3
    segments, num_pix_pkg = build_histogram(vec(intensity_array), num_segments)    
    # remove the first elements 
    # counts = counts[1:end]
    # compute segments 
    segments = compute_histogram(dcm, collect(segments)[2:end])

    # get segments indexes 
    seg_indexs = getindexes.(segments) 
    # get num pixels inside each segment
    num_pixes = getnumpix.(segments)
    # test results
    #TODO: fix num_pix_pkg == num_pixes
    @test num_pix_pkg[end] == num_pixes[end]
    # get segments intensity vector 
    int_segments = getint.(segments)
    # compute means for each segments  
    mean_segments = mean.(segments)
    # compute std for each segments  
    std_segments = std.(segments)
end