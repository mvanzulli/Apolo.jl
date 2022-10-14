#################################################
# Main types and functions to deal with images  #
#################################################

module Images

# Import external dependencies to overload
import Statistics: mean, std
import Base: size, length
# Import local dependencies to overload
import ..ForwardProblem: getgrid, getdim

# Add libraries to use
using DICOM: DICOMData, dcmdir_parse, @tag_str
using SparseArrays: spzeros

# Export interface functions and types
export AbstractImage, MedicalImage, MedicalSegment,
    load_dicom, getintensity, getnumpix,
    getspacing, gethyp, getorientation, 
    get_patient_name, getorigin, build_coordinates,
    getindexes, getnumpix, mean,
    std, getvarcoef, compute_histogram


# Scalar parameters for the interface
const ORIENTATION_TOLERANCE = 1e-2 # this parameter is the tolerance to not consider RAS to IJK rotation

#TODO: med_img(2,2,4) en pixeles
# med_img[1,2,3] = 

""" Abstract supertype for all images.

The following methods are provided by the interface:

- `getintensity(img)`      -- returns the intensity array of the image `img`. 
- `get_intensity_type(img)`-- returns the intensity data type.
- `getnumpix(img)`         -- returns the image resolution in pixels. 
- `getspacing(img)`        -- returns the space between pixels. 
- `getorigin(img)`         -- returns the image origin coordinates. 
- `getlength(img)`         -- returns the image dimensions in metric units.
"""
abstract type AbstractImage{T,D} end 

## Methods for abstract parameters:
" Gets the image intensity array"
getintensity(img::AbstractImage) = img.intensity

" Gets the image dimension"
getdim(::AbstractImage{T,D}) where {T,D} = D

" Gets the image intensity array data type"
get_intensity_type(::AbstractImage{T,D}) where {T,D} = T

" Gets the image size"
size(img::AbstractImage) = size(getintensity(img))

" Gets the image dimensions in voxels"
getnumpix(img::AbstractImage) = img.dimension

" Gets the image spacing between pixels"
getspacing(img::AbstractImage) = img.spacing

" Gets the image origin"
getorigin(img::AbstractImage) = img.origin

" Gets the image metric dimensions"
getlength(img::AbstractImage) = size(img) .* collect(getspacing(img)) 


""" Medical image struct.

### Fields:

- `intensity`   -- intensity array
- `dimension`   -- number of pixels in each direction
- `spacing`     -- space in each direction between each pixel
- `origin`      -- origin coordinates (considering the origin located at [1,1,1] )
- `orientation` -- patient orientation tolerance ORIENTATION_TOLERANCE 
- `hyp_params`  -- last slice DICOM image with its corresponding hyper parameters 

"""
struct MedicalImage{T,D,F} <:AbstractImage{T,D}
    intensity::Array{T,D}
    dimension::NamedTuple
    spacing::NamedTuple
    origin::NTuple{D,<:Real}
    orientation::Symbol
    hyp_params
end
#AbstractImage{D}
#TODO: Fake Image 
#TODO: VTKImage
#TODO:Cashear al medical image PointEvalHandler  
#TODO implement interpolation 
# (med_img::MedicalImage)(x, y, z) = # llamada al eval handler 

" Gets the image hyper-parameters"
getdim(med_img::MedicalImage) = med_img.dimension
" Gets the image hyper-parameters"
gethyp(med_img::MedicalImage) = med_img.hyp_params


" Gets the image coordinates"
function build_coordinates(med_img::MedicalImage{T,2}) where T
    
    # get origin
    (oₓ, oⱼ, oₖ)  = getorigin(med_img)
    # get spacing
    (Δₓ, Δⱼ, Δₖ)  = getspacing(med_img) |> collect
    # get length 
    (Lₓ, Lⱼ, Lₖ) = getlength(med_img)
    # get length 
    (nₓ, nⱼ, nₖ) = getnumpix(med_img) |> collect
    
    # build start and finish 
    x_begin = oₓ + Δₓ/2; x_end =  Lₓ - Δₓ/2
    y_begin = oⱼ + Δⱼ/2; y_end =  Lⱼ - Δⱼ/2
    z_begin = oₖ + Δₖ/2; z_end =  Lₖ - Δₖ/2
    
    # build coordinates 
    x_coords = LinRange(x_begin, x_end, nₓ)
    y_coords = LinRange(y_begin, y_end, nⱼ)
    z_coords = LinRange(z_begin, z_end, nₖ)

    return x_coords, y_coords, z_coords
end
"Gets the patient name"
get_patient_name(med_img::MedicalImage) = gethyp(med_img)[tag"PatientName"]

"Extracts the pixel by its index"#TODO
# Base.getindex(m_img::MedicalImage, i1::Int, i2::Int, i3::Int) = m_img.intensity[i1, i2, i3]


" Extract the pixel data form a DICOM array"
function getpixeldata(dcm_array)
    if length(dcm_array) == 1
        return only(dcm_array).PixelData
    else
        return cat([dcm.PixelData for dcm in dcm_array]...; dims=3)
    end
end

function getorientation(dcm::MedicalImage{T,D,DICOMData}) where {T,D}
    # extract DICOM hyper parameters
    hyp = gethyp(dcm)
    # run 
    orient = getorientation(hyp)
    return orient
end
"Extract the image orientation from a DICOMData"
function getorientation(dcm::DICOMData)
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
getspacing(dcm) = (sagital=dcm[tag"PixelSpacing"][1], coronal=dcm[tag"PixelSpacing"][2], axial=dcm[tag"SliceThickness"])

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
    intensity = getpixeldata(dcms)
    # extract the intensity array data type
    intensity_type = eltype(intensity)
    # extract the image spacing 
    spacing = getspacing(dcm)
    # extract the image dimension 
    dimension = (sagital=size(intensity, 1), coronal=size(intensity, 2), axial=size(intensity, 3))
    dim_number = length(dimension)
    # extract the image orientation
    orientation = getorientation(dcm)
    # construct and build the image 
    med_img = MedicalImage{intensity_type, dim_number, format}(intensity, dimension, spacing, (0, 0, 0), orientation, dcm)
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
getindexes(med_seg::MedicalSegment) = med_seg.indexes
" Gets the number of pixels inside an image segment"
getnumpix(med_seg::MedicalSegment) = med_seg.num_pix
" Gets the image segment intensity array"
getintensity(med_seg::MedicalSegment) = med_seg.intensity_vec

"Computes the histdogram given a vector of intensity edges"
function compute_histogram(
    med_img::MedicalImage,
    edges::Vector{<:Number},
    ) 
    # get the intensity array
    intensity = getintensity(med_img)

    # sort edges min to max value
    edges = sort(edges)

    # check edges are inside the image range 
    (edgeₘᵢₙ, edgeₘₐₓ) = extrema(edges) 
    (intensityₘᵢₙ, intensityₘₐₓ) = extrema(intensity) 
    inside_bool = edgeₘᵢₙ ≥ intensityₘᵢₙ &&  edgeₘₐₓ ≤ intensityₘₐₓ 
    !inside_bool && throw(BoundsError("intensity edges must be inside img_int"))

    # number of  edges
    num_segments = length(edges) + 1

    # initialize
    segments = Vector{MedicalSegment}(undef, num_segments )

    # compute the number of voxels that are inside the borders
    for segment_idx in 1:num_segments
        # segment borders
        # left
        borderₘᵢₙ = if segment_idx == 1
            intensityₘᵢₙ
        else      
            edges[segment_idx - 1]
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
mean(med_seg::MedicalSegment) = mean(getintensity(med_seg))
" Gets the intensity standard deviation of an image segment"
std(med_seg::MedicalSegment) = std(getintensity(med_seg))
" Gets the variation coefficient of an image segment"
varcoef(med_seg::MedicalSegment) = getstd(med_seg) / getmean(med_seg)


end # end module

