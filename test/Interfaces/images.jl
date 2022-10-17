################
# Images tests #
################

# Using internal packages to test 
using Apolo:create_grid,dimension
using Apolo.Images


# Using external packages to test 
using Test: @test, @testset
using AutoHashEquals
using HistogramThresholding: build_histogram
using Statistics:  var, mean
using DICOM: @tag_str


@testset "Generic images" begin
    
    # 2D image 
    # image dimensions 
    spacing_img = (.5, .25)
    num_pixels = (4, 3)
    start_img = (1.,1.)
    image_dimension = length(spacing_img)
    finish_img = start_img .+ num_pixels .* spacing_img
    length_img = finish_img .- start_img
    # image intensity
    intensity_vec = collect(1:12)
    intensity_mat = reshape(intensity_vec, num_pixels )
    # create the image grid
    grid_img = create_img_grid(start_img, spacing_img, length_img, num_pixels)

    # generate image 
    generic_img = GenericImage(intensity_mat, num_pixels, spacing_img, start_img)

    # test extractors
    @test intensity(generic_img) == intensity_mat
    @test intensity_type(generic_img) == eltype(intensity_mat)
    @test dimension(generic_img) == image_dimension
    @test numpix(generic_img) == num_pixels
    @test spacing(generic_img) == spacing_img
    @test start(generic_img) == start_img
    @test finish(generic_img) == finish_img
    @test length(generic_img) == length_img

    # test coordinates generator
    coords = coordinates(generic_img)
    @test spacing_img[1] ==  coords[1][2] - coords[1][1] 
    @test spacing_img[2] ==  coords[2][2] - coords[2][1] 
    extremaᵢ, extremaⱼ = extrema.(coords)
    @test extremaᵢ == (start_img[1] + spacing_img[1]/2, finish_img[1] -  spacing_img[1]/2)
    @test extremaⱼ == (start_img[2] + spacing_img[2]/2, finish_img[2] -  spacing_img[2]/2)

    # test mesh generator
    grid_img_extracted = getgrid(generic_img)
    @test grid_img.cells == grid_img_extracted.cells
    @test grid_img.nodes == grid_img_extracted.nodes
    @test grid_img.facesets == grid_img_extracted.facesets

end



# Constant variables and tpyes
# this file is executed from /test/
const DCM_TO_LOAD = "./DICOMImages/"


# this tests is also validated with 3D slicer results
@testset "Reading DICOM" begin
    # load dicom images
    dcm = load_dicom(DCM_TO_LOAD)

    # get image dimensions
    dim = getdim(dcm);
    num_pix = (256, 256, 400)
    named_num_pix = (sagital = num_pix[1], coronal = num_pix[2], axial = num_pix[3])
    @test dim == named_num_pix

    # get image intensity array
    intensity = getintensity(dcm)
    @test size(intensity) == num_pix

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


