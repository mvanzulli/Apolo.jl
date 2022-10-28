##########################################################
# Main types and functions to handle with medical images #
##########################################################


import Statistics: mean, std
import ..Images: finish_grid, intensity, intensity_type, num_pix, start_grid, spacing

using DICOM: DICOMData, dcmdir_parse, @tag_str
using Apolo: create_rectangular_grid, ⊂

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
hyper_parameters(med_img::MedicalImage) = med_img.hyp_params


"Gets the patient name"
patient_name(med_img::MedicalImage) = hyper_parameters(med_img)[tag"PatientName"]

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
    hyp = hyper_parameters(dcm)
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
spacing(dcm::DICOMData) = (sagital=dcm[tag"PixelSpacing"][1], coronal=dcm[tag"PixelSpacing"][2], axial=dcm[tag"SliceThickness"])

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
    intensity_dcm = pixeldata(dcms)
    # extract the intensity array data type
    intensity_type_dcm = eltype(intensity_dcm)
    # extract the image spacing
    spacing_dcm = spacing(dcm)
    # extract the image dimension
    dimension_dcm = (
        sagital=size(intensity_dcm, 1),
        coronal=size(intensity_dcm, 2),
        axial=size(intensity_dcm, 3)
    )

    dim_number = length(dimension_dcm)
    # extract the image orientation
    orientation_dcm = orientation(dcm)
    # construct and build the image
    med_img = MedicalImage{intensity_type_dcm,dim_number,format}(
        intensity_dcm,
        dimension_dcm,
        spacing_dcm,
        (0, 0, 0),
        orientation_dcm,
        dcm)
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
numpix(med_seg::MedicalSegment) = med_seg.num_pix
" Gets the image segment intensity array"
intensity(med_seg::MedicalSegment) = med_seg.intensity_vec

"Computes the histdogram given a vector of intensity edges"
function compute_histogram(
    med_img::MedicalImage,
    edges::Vector{<:Number},
)
    # get the intensity array
    intensity_med = intensity(med_img)

    # sort edges min to max value
    edges = sort(edges)

    # check edges are inside the image range
    (edgeₘᵢₙ, edgeₘₐₓ) = extrema(edges)
    (intensityₘᵢₙ, intensityₘₐₓ) = extrema(intensity_med)
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
        indexes_segment = findall(is_in_segment, intensity_med)

        # compute the intensity segment array
        intensity_segment = intensity_med[indexes_segment]

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
