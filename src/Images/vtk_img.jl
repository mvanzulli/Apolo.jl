#######################################################
# Main types and functions to handle with .VTK images #
#######################################################

import ..Images: finish_grid, intensity, intensity_type, num_pixels, start_grid, spacing
import ..Geometry:  ⊂, ⊄, coordinates, cartesian_index, dimension, extrema, finish,
grid, start
import ..Geometry: _interpolate, _extrapolate

using ..Geometry: AbstractStructuredGrid
using ..Geometry: _create_ferrite_rectangular_grid
using ..Images: AbstractIntensity
using DICOM: DICOMData, dcmdir_parse, @tag_str

export VTKImage, path

""" VTK image struct.

### Fields:

- `intensity`   -- image intensity
- `dimension`   -- number of pixels in each direction
- `spacing`     -- space in each direction between each pixel
- `start`      -- start coordinates (considering the start located at [1,1,1] )

"""
struct VTKImage{D,T,I<:AbstractIntensity,G<:AbstractStructuredGrid} <: AbstractImage{D,T,G}
    intensity::I
    num_pixels::NTuple{D,<:Int}
    start::NTuple{D,T}
    spacing::NTuple{D,T}
    grid::G
    path::String
end

"VTKImage constructor"
function VTKImage(
    intensity_array::Array{T,D},
    spacing_img::NTuple{D,T} = Tuple(ones(T,D)),
    start_img::NTuple{D,T} = Tuple(zeros(T,D)),
    path_img::String = "";
    ferrite_grid::Bool = true,
) where {T,D}

    ferrite_grid == false && throw(ArgumentError("The grid must be a ferrite subtype"))

    # number of pixels
    num_pixels = size(intensity_array)

    # compute end point
    finish_img = finish(start_img, num_pixels, spacing_img)
    length_img = length(start_img, finish_img)

    # create grid
    fgrid = create_ferrite_img_fgrid(start_img, spacing_img, length_img, num_pixels)

    # convert intensity to ferrite nomenclature
    fintensity= FerriteIntensity(intensity_array, fgrid)

    # instantiate generic grid
    VTKImage(fintensity, num_pixels, start_img, spacing_img, fgrid, path_img)

end

path(vtk_img::VTKImage) = vtk_img.path_img
