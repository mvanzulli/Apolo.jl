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
using Apolo: create_rectangular_grid, ⊂
using Ferrite: Vec, Grid, Quadrilateral, CellIterator,
    DofHandler, PointEvalHandler, Lagrange, RefCube,
    get_point_values, close!, getnodes

# Export interface functions and types
export AbstractImage, MedicalImage, MedicalSegment, GenericImage,
    # Extractors AbstractImage
    intensity, intensity_type, numpix, total_numpix,
    spacing, start, finish, length, coordinates,
    # Features AbstractImage 
    create_img_grid,
    # Extractors MedicalImage
    hyper_parameters, orientation, patient_name, indexes,
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

" Gets the total number of pixels"
total_numpix(img::AbstractImage) = prod(numpix(img))

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


"Checks if a index is inbounds the image dimensions"
_index_is_inbounds(onedim_index::Int, img::AbstractImage) = 1 ≤ onedim_index ≤ total_numpix(img)

"Checks if a index is inbounds the image dimensions"
function _index_is_inbounds(ndim_index::CartesianIndex{D}, img::AbstractImage{T,D}) where {T,D} 
    
    num_pixels = numpix(img)
    
    # Check all the indexes are inbounds
    index_is_in = true 
    for direction in 1:D 
        index_is_in = index_is_in && 1 ≤ ndim_index[direction] ≤ num_pixels[direction]
    end

    return index_is_in

end

"Checks if a point is inbounds the image "
function _point_is_inbounds(p::NTuple{D}, img::AbstractImage{T,D}) where {T,D} 
    
    # Extract borders
    start_point = start(img)
    finish_point = finish(img)

    p_is_in = true

    for axis in 1:D
        p_is_in = p_is_in &&  start_point[axis] ≤ p[axis] ≤ finish_point[axis]
    end

    return p_is_in
end


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
    num_elements = num_pixels .- Tuple(ones(Int, D))
    # compute start point
    start_mesh = start_img .+ spacing_img ./ 2
    # compte finish point
    finish_mesh = finish_img .- spacing_img ./ 2

    # generate grid
    grid = create_rectangular_grid(D, num_elements, start_mesh, finish_mesh, Quadrilateral)

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

(img::AbstractImage{T,2} where {T})(x, y) = _eval_intensity(img, (x, y))
(img::AbstractImage{T,3} where {T})(x, y, z) = # llamada al eval handler 

"Internal function to evaluate the image intensity at a generic point via Ferrite Grid"
function _eval_intensity(img::AbstractImage{T,2}, p::NTuple{2,<:Real}) where {T}

    # Check p is inisde the image
    _point_is_inbounds(p, img)
    # Check if p is inside the grid image 
    if _is_inside_img_grid(img, p) 
        intensity_p = _eval_intensity_inside_grid(img, p)
    else
        # Extrapolate if is not inside the grid 
    end
    return intensity_p
end

"Dof Handler constructor given an image "
function DofHandler(img::AbstractImage)
    
    # create dof handler
    grid_img = getgrid(img)
    dh_img = DofHandler(grid_img)

    # push intensity field
    push!(dh_img, :intensity, 1, Lagrange{2,RefCube,1}())
    close!(dh_img)

    return dh_img

end


"Internal function to evaluate the image at a generic point (inside the image grid) "
function _eval_intensity_inside_grid(
    img::AbstractImage{T,2},
    p::NTuple{2,<:Real}
) where {T}
    #TODO: Compute the ferrite intensity array only one time not every time this function is called
    
    # get the image intensity and transform it into a vector with ferrite.jl nomenclature
    grid_img = getgrid(img)
    intensity_ferrite_vec = _build_ferrite_img_intensity(img, grid_img)    # create dof handler

    # dof handler
    dh_img = DofHandler(img)

    # create evaluation handlar
    eval_point = [Vec(p)]
    ph_img = PointEvalHandler(grid_img, eval_point)

    # evaluate 
    i_points = get_point_values(
        ph_img, dh_img,
        intensity_ferrite_vec,
        :intensity
    )

    return i_points
end

" Build intensity vector according to ferrite nomenclature (cell by cell CCW)"
function _build_ferrite_img_intensity(img::AbstractImage, grid::Grid)
    
    # ferrite intensity vector 
    intensity_ferrite_vec = Vector{Float64}()

    # extract image intensity array
    intensity_img = intensity(img)

    # since cells share nodes (nodes that has been added are redundant)
    added_nodes_index = Vector{Int}() 

    # iterate over each grid cell
    for cell in CellIterator(grid)
        
        # get cell nodes
        nodes_cell = getnodes(cell)

        # add intensity if the node is not been added yet
        for node_cell in nodes_cell

            # check if is already been added
            if node_cell ∉ added_nodes_index  
                push!(added_nodes_index, node_cell)
                push!(intensity_ferrite_vec, intensity_img[_node_cartesian_index(node_cell, img)])
            end
        
        end

    end

    return intensity_ferrite_vec
end

"Computes the Cartesian index given a 1D index node. "
function _node_cartesian_index(onedim_index::Int, img::AbstractImage)
    # Node indexes in each direction respectively
    num_pixels = collect(numpix(img))

    # Check index_nodes is inbounds 
    !(_index_is_inbounds(onedim_index, img)) &&
     throw(ArgumentError(BoundsError("index_node is not inside the image size $(size(img))")))

    # Compute indexes depending on the image dimension 
    if length(num_pixels) == 2 # 2D image
    
        # if rem(onedim_index, num_pix_x) == 0 then index_x = num_pix_x
        (index_x, index_y) = _my_div_rem(onedim_index, num_pixels[1], num_pixels[1])

        return CartesianIndex(index_x, index_y)
        
    elseif length(num_pixels) == 3 # 3D image

        # if rem(onedim_index,num_pix_x * num_pix_y) == 0 then index_z = num_pixels[3]
        (index_z, rem_xy) = _my_div_rem(onedim_index,  num_pixels[1] * num_pixels[2], num_pixels[3])

        # if rem(onedim_index, rem_xy) == 0 then index_x = num_pix_x
        (index_x, index_y) = _my_div_rem(rem_xy, num_pixels[1], num_pixels[1])

        # return the index
        return CartesianIndex(index_x, index_y, index_z)
    
    end

end

function _my_div_rem(onedim_index::Int, batch::Int, edge_val::Int)
    # div rem
    index, rem = divrem(onedim_index, batch)
    
    # case onedim_index ≤ batch
    if index == 0
        return (onedim_index, 1)
    end
    # fix edge case 
    if rem == 0
        index = index 
        rem = edge_val 
        return rem, index
    else
        return rem, index + 1
    end

end

"Checks if p is inside the image grid. "
function _is_inside_img_grid(img, p::NTuple{T,<:Real}) where {T}
    # extract grid 
    grid_img = getgrid(img)
    # check and return
    return p ⊂ grid_img
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


end # end module

