##########################################################
# Main types and functions to handle with ferrite images #
##########################################################

using Ferrite: Quadrilateral, Hexahedron

using ..Geometry: FerriteStructuredGrid
using ..Geometry: _convert_to_ferrite_nomenclature
using ..Images: finish_grid, finish, start_grid, start
using ..Images: _grid_num_elements

import ..Images: intensity, value

const FSG = FerriteStructuredGrid

export FerriteImage, FerriteIntensity, create_ferrite_img_fgrid


""" Ferrite intensity struct.

### Fields:

- `intensity`   -- Intensity array with ferrite nomenclature
- `indexes`-- Vector of cartesian indexes that maps ferrite intensity value
"""
struct FerriteIntensity{D,T} <: AbstractIntensity{D,T}
    intensity_vec::Vector{T}
    indexes::Vector{CartesianIndex{D}}
end

_intensity_vec(fint::FerriteIntensity) = fint.intensity_vec
_indexes(fint::FerriteIntensity) = fint.indexes

function FerriteIntensity(
    intensity_array::Array{T,D},
    fgrid::FerriteStructuredGrid{D}
) where {T,D}

    # convert intensity to ferrite nomenclature
    int_ferrite, indexes = _convert_to_ferrite_nomenclature(intensity_array, fgrid)

    return FerriteIntensity{D,T}(int_ferrite, indexes)
end

""" Ferrite image struct.

### Fields:

- `intensity`  -- Intensity array with ferrite nomenclature
- `num_pixels` -- Number of pixels in each direction
- `start`      -- Start coordinates (considering the start located at [1,1,1] )
- `spacing`    -- Space in each direction between each pixel
- `grid`       -- Ferrite grid inside the image
"""
struct FerriteImage{D,T} <: AbstractImage{D,T}
    fintensity::FerriteIntensity{D,T}
    num_pixels::NTuple{D,<:Integer}
    start::NTuple{D,<:Real}
    spacing::NTuple{D,<:Real}
    grid::FSG
end

"Creates a FerriteStructuredGrid inside a 2D image frame"
function create_ferrite_img_fgrid(
    start_img::NTuple{2,<:Real},
    spacing_img::NTuple{2,<:Real},
    length_img::NTuple{2,<:Real},
    num_pixels::NTuple{2,<:Integer},
)

    # Image parameters
    finish_img = start_img .+ length_img

    # Mesh parameters
    num_elements = _grid_num_elements(num_pixels)
    start_grid_point = start_grid(start_img, spacing_img)
    finish_grid_point = finish_grid(finish_img, spacing_img)
    element = Quadrilateral

    # generate grid
    fgrid = FerriteStructuredGrid(start_grid_point, finish_grid_point, num_elements, element)

    return fgrid
end

"Creates a FerriteStructuredGrid inside a 3D image frame"
function create_ferrite_img_fgrid(
    start_img::NTuple{3,<:Real},
    spacing_img::NTuple{3,<:Real},
    length_img::NTuple{3,<:Real},
    num_pixels::NTuple{3,<:Integer},
)

    # Image parameters
    finish_img = start_img .+ length_img

    # Mesh parameters
    num_elements = _grid_num_elements(num_pixels)
    start_grid_point = start_grid(start_img, spacing_img)
    finish_grid_point = finish_grid(finish_img, spacing_img)
    element = Hexahedron

    # generate grid
    fgrid = FerriteStructuredGrid(start_grid_point, finish_grid_point, num_elements, element)

    return fgrid
end

# Cash the Ferrite.Grid into grid field of MedicalImage
function FerriteImage(
    intensity_array::Array{T,D},
    num_pixels::NTuple{D,<:Integer},
    start_img::NTuple{D,<:Real},
    spacing_img::NTuple{D,<:Real},
) where {T,D}

    # compute end point
    finish_img = finish(start_img, num_pixels, spacing_img)
    length_img = length(start_img, finish_img)

    # create grid
    fgrid = create_ferrite_img_fgrid(start_img, spacing_img, length_img, num_pixels)

    # convert intensity to ferrite nomenclature
    fintensity= FerriteIntensity(intensity_array, fgrid)

    # instantiate generic grid
    return FerriteImage(fintensity, num_pixels, start_img, spacing_img, fgrid)
end

function intensity(fimg::FerriteImage{D,T}) where {D,T}

    fintensity = fimg.fintensity
    intensity_vec = _intensity_vec(fintensity)
    indexes = _indexes(fintensity)

    intensity_array = Array{T,D}(undef, num_pixels(fimg))

    [intensity_array[indexes[i]] = int for (i, int) in enumerate(intensity_vec)]

    return intensity_array
end
