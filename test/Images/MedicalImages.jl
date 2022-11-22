#######################
# Medical Image Tests #
#######################
using Test: @testset, @test
using Apolo.Geometry.FerriteGrids: _convert_to_ferrite_nomenclature
using Apolo.Images
using Apolo.Images.MedicalImages
using Apolo.Images: _total_num_pixels, _index_is_inbounds, _eval_intensity,
    _interpolate, _extrapolate

const INTERVAL_START = LinRange(-10.0, 10.0, 20)
const INTERVAL_POS = LinRange(1.0, 10.0, 20)
const INTERVAL_OFFSET = LinRange(-10.0, -1.0, 20)
const TOLERANCE = 1e-3
const DCM_TO_LOAD = "./test/DICOMImages/"


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
