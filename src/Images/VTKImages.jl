"""
Module defining vtk images (an image implemented reed from a .VTK file).
"""
module VTKImages

using Apolo.Geometry: AbstractStructuredGrid
using Apolo.Geometry.FerriteGrids: FerriteStructuredGrid
using ..Images: AbstractImage, AbstractIntensity, FerriteIntensity
using ..Images: create_ferrite_img_fgrid, finish

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

"VTKImage constructor given an `intensity_array`, an spacing vector between pixels `spacing_img`,
start point `start_img`, vtk image path `path_img` and a given `ferrite_grid` boolean ."
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


"VTKImage constructor given an `intensity_array`, an spacing vector between pixels `spacing_img`,
start point `start_img`, vtk image path `path_img` and a given ferrite grid `fgrid`."
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

    # TODO: Check grid dimensions matches start and spacing

    # instantiate generic grid
    return VTKImage(fintensity, num_pixels, start_img, spacing_img, fgrid, path_img)

end

end #endmodule
