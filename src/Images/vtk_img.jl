#######################################################
# Main types and functions to handle with .VTK images #
#######################################################

using ..Geometry: AbstractStructuredGrid, FerriteStructuredGrid
using ..Images: AbstractImage, AbstractIntensity
using ..Images: create_ferrite_img_fgrid

export VTKImage

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

"VTKImage constructor given an intensity array"
function VTKImage(
    intensity_array::Array{T,D},
    spacing_img::NTuple{D,T}=Tuple(ones(T, D)),
    start_img::NTuple{D,T}=Tuple(zeros(T, D)),
    path_img::String="";
    ferrite_grid::Bool=true
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
    fintensity = FerriteIntensity(intensity_array, fgrid)

    # instantiate generic grid
    return VTKImage(fintensity, num_pixels, start_img, spacing_img, fgrid, path_img)

end

function VTKImage(
    intensity_array::Array{T,D},
    fgrid::FerriteStructuredGrid,
    spacing_img::NTuple{D,T}=Tuple(ones(T, D)),
    start_img::NTuple{D,T}=Tuple(zeros(T, D)),
    path_img::String="";
) where {T,D}

    # number of pixels
    num_pixels = size(intensity_array)

    # convert intensity to ferrite nomenclature
    fintensity = FerriteIntensity(intensity_array, fgrid)

    # instantiate generic grid
    return VTKImage(fintensity, num_pixels, start_img, spacing_img, fgrid, path_img)

end
