#################################################
# Main types and functions to deal with images  #
#################################################

module Images

# Import external dependencies to overload
import Statistics: mean, std
# Import local dependencies to overload
import ..ForwardProblem: getdim, getgrid

# Add libraries to use
using DICOM: dcmdir_parse, @tag_str
using SparseArrays: spzeros

# Export interface functions and types
export MedicalImage, MedicalSegment,
    load_dicom, getintensity, getdim,
    getspacing, gethyp, getorientation,
    getorigin, 
    getindexes, getnumpix, getint, mean,
    std, getvarcoef, compute_histogram


# Scalar parameters for the interface
const ORIENTATION_TOLERANCE = 1e-2 # this parameter is the tolerance to not consider RAS to IJK rotation

""" Medical image struct.

### Fields:

- `intensity`   -- intensity array
- `dimension`   -- number of pixels in each direction
- `spacing`     -- space in each direction between each pixel
- `origin`      -- origin coordinates (considering the origin located at [1,1,1] )
- `orientation` -- patient orientation tolerance ORIENTATION_TOLERANCE 
- `hyp_params`  -- last slice DICOM image with its corresponding hyper parameters 

"""
struct MedicalImage{T}
    intensity::Array{<:Number,3}
    dimension::NamedTuple
    spacing::NamedTuple
    origin::NTuple{3,<:Real}
    orientation::Symbol
    hyp_params
end

" Gets the image intensity array"
getintensity(med_img::MedicalImage) = med_img.intensity
" Gets the image dimensions"
getdim(med_img::MedicalImage) = med_img.dimension
" Gets the image spacing between pixels"
getspacing(med_img::MedicalImage) = med_img.spacing
" Gets the image orientation"
getorientation(med_img::MedicalImage) = med_img.orientation
" Gets the image origin"
getorigin(med_img::MedicalImage) = med_img.origin
" Gets the image hyper-paramters"
gethyp(med_img::MedicalImage) = med_img.hyp_params

"Extracts the pixel by its index"#TODO
# Base.getindex(m_img::MedicalImage, i1::Int, i2::Int, i3::Int) = m_img.intensity[i1, i2, i3]


" Extract the pixel data form a DICOM  array"
function getpixeldata(dcm_array)
    if length(dcm_array) == 1
        return only(dcm_array).PixelData
    else
        return cat([dcm.PixelData for dcm in dcm_array]...; dims=3)
    end
end

"Extract the image orientation"
function getorientation(dcm)
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
function load_dicom(dir)
    # load the dicom directory
    dcms = dcmdir_parse(dir)
    # 'dcms' only one type of serie is admitted so we have to check that 
    unique_serie = unique([dcm.SeriesInstanceUID for dcm in dcms])
    length(unique_serie) > 1 && throw(ArgumentError("The SeriesInstanceUID of each image in $dir must be the same"))
    dcm = dcms[1]
    #  Construct the named medical image
    # extract the image type
    img_type = typeof(dcm)
    # extract the intensity array
    intensity = getpixeldata(dcms)
    # extract the image spacing 
    spacing = getspacing(dcm)
    # extract the image dimension 
    dimension = (sagital=size(intensity, 1), coronal=size(intensity, 2), axial=size(intensity, 3))
    # extract the image orientation
    orientation = getorientation(dcm)
    # construct and build the image 
    return MedicalImage{img_type}(intensity, dimension, spacing, (0, 0, 0), orientation, dcm)
end


""" Medical segment image struct. This corresponds to a segment of intensity in 
the histogram sorted by intensity values. 

### Fields:

- `indexes`   -- indexes of the image in i,j,k coordinates that corresponds to this segment.
- `num_pix`   -- number of pixels inside this segment.
- `intensity` -- segment intensity array.
"""
struct MedicalSegment
    indexes::Vector{CartesianIndex}
    num_pix::Integer
    intensity_vec::Array{<:Number}
end


" Gets medical image indexes"
getindexes(med_seg::MedicalSegment) = med_seg.indexes
" Gets the number of pixels inside an image segment"
getnumpix(med_seg::MedicalSegment) = med_seg.num_pix
" Gets the image segment intensity array"
getint(med_seg::MedicalSegment) = med_seg.intensity_vec

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
    (int_imgₘᵢₙ, int_imgₘₐₓ) = extrema(intensity) 
    inside_bool = edgeₘᵢₙ ≥ int_imgₘᵢₙ &&  edgeₘₐₓ ≤ int_imgₘₐₓ 
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
            int_imgₘᵢₙ
        else      
            edges[segment_idx - 1]
        end 
        #right
        borderₘₐₓ = if segment_idx == num_segments
            int_imgₘₐₓ
        else
            edges[segment_idx]
        end
        # anonymmus function to check if it inside
        is_in_segment(x) = borderₘᵢₙ ≤ x ≤ borderₘₐₓ

        # compute the cartesian indexes of pixels that satisfy is_in_segment condition  
        indexes_segment = findall(is_in_segment, intensity)
        
        # compute the intensity segment array 
        int_segment = intensity[indexes_segment]

        # number of pixels inside the segment
        counts_segment = length(indexes_segment)

       # build segment 
        segment_obj = 
        MedicalSegment(
            indexes_segment,
            counts_segment,
            int_segment,
            )

        # fill segments vector
        segments[segment_idx] = segment_obj
    end

    return segments
end


" Gets the image segment intensity mean"
mean(med_seg::MedicalSegment) = mean(getint(med_seg))
" Gets the intensity standard deviation of an image segment"
std(med_seg::MedicalSegment) = std(getint(med_seg))
" Gets the variation coefficient of an image segment"
varcoef(med_seg::MedicalSegment) = getstd(med_seg) / getmean(med_seg)


end # end module

