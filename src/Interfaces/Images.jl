#################################################
# Main types and functions to deal with images  #
#################################################

module Images

# Import external dependencies to overload
import Statistics: mean, std
import Base: size, length

# Import local dependencies to overload
import ..ForwardProblem: dimension, getgrid

# Add libraries to use
using AutoHashEquals
using DICOM: DICOMData, dcmdir_parse, @tag_str
using Ferrite: Grid, Quadrilateral, DofHandler, get_point_values
using Apolo: create_grid

# Export interface functions and types
export AbstractImage, MedicalImage, MedicalSegment, GenericImage,
    # Extractors AbstractImage
    intensity, intensity_type, numpix,
    spacing, start, finish, length, coordinates,
    # Features AbstractImage 
    create_img_grid,
    # Extractors MedicalImage
    hyper_params, orientation, patient_name, indexes,
    # Medical image features
    mean, load_dicom, std, getvarcoef, compute_histogram


# Scalar parameters for the interface
const ORIENTATION_TOLERANCE = 1e-2 # this parameter is the tolerance to not consider RAS to IJK rotation

""" Abstract supertype for all images.

The following methods are provided by the interface:

- `intensity(img)`         -- returns the intensity array of the image `img`. 
- `intensity_type(img)`    -- returns the intensity data type.
- `numpix(img)`            -- returns the image resolution in pixels. 
- `spacing(img)`           -- returns the space between pixels. 
- `start(img)`             -- returns the start point coordinates. 
- `finish(img)`            -- returns the finish point coordinates.
- `length(img)`            -- returns the image dimensions in metric units.
- `coordinates(img)`       -- returns the image coordinates in three linear ranges.
- `grid(img)`              -- returns the image grid.
"""
abstract type AbstractImage{T,D} end

## Methods for abstract parameters:
" Gets the image intensity array"
intensity(img::AbstractImage) = img.intensity

" Gets the image dimension"
dimension(::AbstractImage{T,D}) where {T,D} = D

" Gets the image intensity array data type"
intensity_type(::AbstractImage{T,D}) where {T,D} = T

" Gets the image size"
size(img::AbstractImage) = size(intensity(img))

" Gets the image dimensions in voxels"
numpix(img::AbstractImage) = img.num_pixels

" Gets the image spacing between pixels"
spacing(img::AbstractImage) = img.spacing

" Gets the image start point"
start(img::AbstractImage) = img.start

" Gets the image metric dimensions"
length(img::AbstractImage) = size(img) .* spacing(img)

" Gets the image finish point"
finish(img::AbstractImage) = start(img) .+ length(img)

" Gets the image grid"
getgrid(img::AbstractImage) = img.grid

" Gets the image coordinates"
function coordinates(img::AbstractImage)

    # get start (start point)
    start_img = start(img)
    # get spacing
    spacing_img = spacing(img) |> collect
    # get length 
    length_img = length(img)
    # get number of pixels 
    num_pixels = numpix(img) |> collect
    # compute finish 
    finish_img = start_img .+ length_img


    # build start and finish point of the pixels grid
    X_begin = start_img .+ spacing_img ./ 2
    X_end = finish_img .- spacing_img ./ 2

    # build coordinates 
    coords = LinRange.(X_begin, X_end, num_pixels)

    return coords
end

"Creates mesh grid inside the image frame"
function create_img_grid(
    start_img::NTuple{D,<:Real},
    spacing_img::NTuple{D,<:Real},
    length_img::NTuple{D,<:Real},
    num_pixels::NTuple{D,<:Integer},
) where {D}

    # Image parameters
    finish_img = start_img .+ length_img

    # Mesh parameters
    num_elements = num_pixels .- Tuple(ones(Int,D))
    # compute start point
    start_mesh = start_img .+ spacing_img ./ 2
    # compte finish point
    finish_mesh = finish_img .- spacing_img ./ 2

    # generate grid
    grid = create_grid(2, num_elements, start_mesh, finish_mesh, Quadrilateral)

    return grid
end


""" Generic image struct.

### Fields:

- `intensity`  -- Intensity array
- `num_pixels` -- Number of pixels in each direction
- `spacing`    -- Space in each direction between each pixel
- `start`      -- Start coordinates (considering the start located at [1,1,1] )
- `grid`       -- Ferrite grid inside the image
"""
struct GenericImage{T,D} <: AbstractImage{T,D}
    intensity::Array{T,D}
    num_pixels::NTuple{D,<:Integer}
    spacing::NTuple{D,<:Real}
    start::NTuple{D,<:Real}
    grid::Grid
    # Cash the Ferrite.Grid into grid field of MedicalImage
    function GenericImage(
        intensity::Array{T,D},
        num_pixels::NTuple{D,<:Integer},
        spacing::NTuple{D,<:Real},
        start::NTuple{D,<:Real},
        ) where {T,D}
        
        # compute end point
        finish = start .+ num_pixels .* spacing
        length_img = finish .- start
        # create grid
        grid = create_img_grid(start, spacing, length_img, num_pixels)
        
        # instantiate generic grid 
        new{T,D}(intensity, num_pixels, spacing, start, grid)
    end
end

(img::AbstractImage{T,2})(x, y, z) = _eval_intenisty(img,(x, y)) 
(img::AbstractImage{T,3})(x, y, z) = # llamada al eval handler 

function _eval_intenisty(img::AbstractImage{T,2}, p::NTuple{2, <:Real}) where T 

    #TODO: Increase performance of this operation
    # get the image intensity and transform it into a vector
    int_ferrite_vec = Vector{Float64}()

    for intᵢ in vec(intensity(img))
        push!(int_ferrite_vec, intᵢ)
    end

    # create dof handler
    grid_img = getgrid(img)
    dh_img = DofHandler(grid_img)

    # push intensity field
    push!(dh_img, :intensity, 1, Lagrange{2,RefCube,1}())
    close!(dh_img)

    # create evaluation handlar
    eval_point = [Vec(p)]
    ph_img = PointEvalHandler(grid_img, eval_point);

    # 
    i_points = Ferrite.get_point_values(
        ph_img, dh_img, 
        int_ferrite_vec,
        :intensity
        )

    return i_points
end


""" Medical image struct.

### Fields:

- `intensity`   -- intensity array
- `dimension`   -- number of pixels in each direction
- `spacing`     -- space in each direction between each pixel
- `start`      -- start coordinates (considering the start located at [1,1,1] )
- `orientation` -- patient orientation tolerance ORIENTATION_TOLERANCE 
- `hyp_params`  -- last slice DICOM image with its corresponding hyper parameters 

"""
struct MedicalImage{T,D,F} <: AbstractImage{T,D}
    intensity::Array{T,D}
    dimension::NamedTuple
    spacing::NamedTuple
    start::NTuple{D,<:Real}
    orientation::Symbol
    hyp_params
    grid::Grid
    # Cash the Ferrite.Grid element of the MedicalImage
    # function MedicalImage{T,D,F}(
    #     intensity::Array{T,D}
    #     dimension::NamedTuple
    #     spacing::NamedTuple
    #     start::NTuple{D,<:Real}
    #     orientation::Symbol
    #     hyp_params) where {T,D,F} 

    # end
end
#AbstractImage{D}
#TODO: Fake Image 
#TODO: VTKImage
#TODO:Cashear al medical image PointEvalHandler  
#TODO implement interpolation 
# (med_img::MedicalImage)(x, y, z) = # llamada al eval handler 

" Gets the image hyper-parameters"
dimension(med_img::MedicalImage) = med_img.dimension
" Gets the image hyper-parameters"
hyper_parameter(med_img::MedicalImage) = med_img.hyp_params


"Gets the patient name"
patient_name(med_img::MedicalImage) = hyper_parameter(med_img)[tag"PatientName"]

"Extracts the pixel by its index"#TODO
# Base.getindex(m_img::MedicalImage, i1::Int, i2::Int, i3::Int) = m_img.intensity[i1, i2, i3]


" Extract the pixel data form a DICOM array"
function pixeldata(dcm_array)
    if length(dcm_array) == 1
        return only(dcm_array).PixelData
    else
        return cat([dcm.PixelData for dcm in dcm_array]...; dims=3)
    end
end

function orientation(dcm::MedicalImage{T,D,DICOMData}) where {T,D}
    # extract DICOM hyper parameters
    hyp = hyper_parameter(dcm)
    # run 
    orient = orientation(hyp)
    return orient
end
"Extract the image orientation from a DICOMData"
function orientation(dcm::DICOMData)
    orientations = Dict(
        :coronal => [1, 0, 0, 0, 0, -1],
        :sagittal => [0, 1, 0, 0, 0, -1],
        :axial => [1, 0, 0, 0, 1, 0]
    )

    # Look for the key whos value is the ImageOrientationPatient field 
    img_orientation_vec = dcm[tag"ImageOrientationPatient"]
    for (orientation, value) in orientations
        if ≈(img_orientation_vec, value, atol=1e-2)
            return orientation
        end
    end
end

"Returns a named tuple with image spacing in axial, sagital and coronal direction "
spacing(dcm) = (sagital=dcm[tag"PixelSpacing"][1], coronal=dcm[tag"PixelSpacing"][2], axial=dcm[tag"SliceThickness"])

" Loads a DICOM directory and returns a named tuple array of different series "
function load_dicom(dir::String)
    # load the dicom directory
    dcms = dcmdir_parse(dir)
    # 'dcms' only one type of serie is admitted so we have to check that 
    unique_serie = unique([dcm.SeriesInstanceUID for dcm in dcms])
    length(unique_serie) > 1 && throw(ArgumentError("The SeriesInstanceUID of each image in $dir must be the same"))
    dcm = dcms[1]
    #  Construct the named medical image
    # extract the image type
    format = typeof(dcm)
    # extract the intensity array
    intensity = pixeldata(dcms)
    # extract the intensity array data type
    intensity_type = eltype(intensity)
    # extract the image spacing 
    spacing = spacing(dcm)
    # extract the image dimension 
    dimension = (sagital=size(intensity, 1), coronal=size(intensity, 2), axial=size(intensity, 3))
    dim_number = length(dimension)
    # extract the image orientation
    orientation = orientation(dcm)
    # construct and build the image 
    med_img = MedicalImage{intensity_type,dim_number,format}(intensity, dimension, spacing, (0, 0, 0), orientation, dcm)
    return med_img
end


""" Medical segment image struct. This corresponds to a segment of intensity in 
the histogram sorted by intensity values. 

### Fields:

- `indexes`   -- indexes of the image in i,j,k coordinates that corresponds to this segment.
- `num_pix`   -- number of pixels inside this segment.
- `intensity` -- segment intensity array.
"""
struct MedicalSegment{T,D}
    indexes::Vector{CartesianIndex{D}}
    num_pix::Integer
    intensity_vec::Array{T}
end


" Gets medical image indexes"
indexes(med_seg::MedicalSegment) = med_seg.indexes
" Gets the number of pixels inside an image segment"
num_pixels(med_seg::MedicalSegment) = med_seg.num_pix
" Gets the image segment intensity array"
intensity(med_seg::MedicalSegment) = med_seg.intensity_vec

"Computes the histdogram given a vector of intensity edges"
function compute_histogram(
    med_img::MedicalImage,
    edges::Vector{<:Number},
)
    # get the intensity array
    intensity = intensity(med_img)

    # sort edges min to max value
    edges = sort(edges)

    # check edges are inside the image range 
    (edgeₘᵢₙ, edgeₘₐₓ) = extrema(edges)
    (intensityₘᵢₙ, intensityₘₐₓ) = extrema(intensity)
    inside_bool = edgeₘᵢₙ ≥ intensityₘᵢₙ && edgeₘₐₓ ≤ intensityₘₐₓ
    !inside_bool && throw(BoundsError("intensity edges must be inside img_int"))

    # number of  edges
    num_segments = length(edges) + 1

    # initialize
    segments = Vector{MedicalSegment}(undef, num_segments)

    # compute the number of voxels that are inside the borders
    for segment_idx in 1:num_segments
        # segment borders
        # left
        borderₘᵢₙ = if segment_idx == 1
            intensityₘᵢₙ
        else
            edges[segment_idx-1]
        end
        #right
        borderₘₐₓ = if segment_idx == num_segments
            intensityₘₐₓ
        else
            edges[segment_idx]
        end
        # anonymmus function to check if it inside
        is_in_segment(x) = borderₘᵢₙ ≤ x ≤ borderₘₐₓ

        # compute the cartesian indexes of pixels that satisfy is_in_segment condition  
        indexes_segment = findall(is_in_segment, intensity)

        # compute the intensity segment array 
        intensity_segment = intensity[indexes_segment]

        # number of pixels inside the segment
        counts_segment = length(indexes_segment)

        # build segment 
        segment_obj =
            MedicalSegment(
                indexes_segment,
                counts_segment,
                intensity_segment,
            )

        # fill segments vector
        segments[segment_idx] = segment_obj
    end

    return segments
end


" Gets the image segment intensity mean"
mean(med_seg::MedicalSegment) = mean(intensity(med_seg))
" Gets the intensity standard deviation of an image segment"
std(med_seg::MedicalSegment) = std(intensity(med_seg))
" Gets the variation coefficient of an image segment"
varcoef(med_seg::MedicalSegment) = getstd(med_seg) / getmean(med_seg)


end # end module

