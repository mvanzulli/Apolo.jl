"""
Module defining ferrte image (an image implemented with a ferrite grid).
"""
module FerriteImages

using Apolo.Geometry.FerriteGrids: FerriteStructuredGrid
using Apolo.Geometry.FerriteGrids: _convert_to_ferrite_nomenclature
using ..Images: AbstractIntensity, AbstractImage
using ..Images: finish_grid, finish, num_pixels, _grid_num_elements, start_grid, start

using Ferrite: Quadrilateral, Hexahedron, DofHandler, PointEvalHandler, Lagrange, RefCube, Vec
using Ferrite: get_point_values, close!

import ..Images: intensity, value

export FerriteImage, FerriteIntensity, create_ferrite_img_fgrid

const FSG = FerriteStructuredGrid


""" Ferrite intensity struct.

### Fields:

- `intensity`   -- Intensity array with ferrite nomenclature.
- `indexes`     -- Vector of cartesian indexes that maps ferrite intensity value.
"""
struct FerriteIntensity{D,T} <: AbstractIntensity{T}
    intensity_vec::Vector{T}
    indexes::Vector{CartesianIndex{D}}
end

"Returns the intensity values of Ferrite Intensity struct `fint`."
value(fint::FerriteIntensity) = fint.intensity_vec

"Returns the indexes of Ferrite Intensity struct `fint`."
_indexes(fint::FerriteIntensity) = fint.indexes

"Returns of a Ferrite Intensity type with an array `intensity_array` and a grid `grid`. "
function FerriteIntensity(
    intensity_array::Array{T,D},
    fgrid::FSG{D}
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
struct FerriteImage{D,T,G<:FSG} <: AbstractImage{D,T,G}
    intensity::FerriteIntensity{D,T}
    num_pixels::NTuple{D,<:Integer}
    start::NTuple{D,<:Real}
    spacing::NTuple{D,<:Real}
    grid::G
end

"Ferrite image constructor given an intensity array `spacing` between pixels and an `start_point`."
function FerriteImage(
    intensity_array::Array{T,D},
    spacing_img::NTuple{D,T}=Tuple(ones(T, D)),
    start_img::NTuple{D,T}=Tuple(zeros(T, D))
) where {T,D}

    num_pixels = size(intensity_array)

    # compute end point
    finish_img = finish(start_img, num_pixels, spacing_img)
    length_img = length(start_img, finish_img)

    # create grid
    fgrid = create_ferrite_img_fgrid(start_img, spacing_img, length_img, num_pixels)

    # convert intensity to ferrite nomenclature
    intensity = FerriteIntensity(intensity_array, fgrid)

    # instantiate generic grid
    FerriteImage{D,T,FSG}(intensity, num_pixels, start_img, spacing_img, fgrid)
end

"Ferrite image constructor given an intensity array `spacing` between pixels and an `start_point`."
function FerriteImage(
    intensity_array::Array{T,D},
    fgrid::G,
    spacing_img::NTuple{D,T}=Tuple(ones(T, D)),
    start_img::NTuple{D,T}=Tuple(zeros(T, D)),
) where {T,D,G<:FSG}

    num_pixels = size(intensity_array)

    # compute end point
    finish_img = finish(start_img, num_pixels, spacing_img)
    length_img = length(start_img, finish_img)

    # convert intensity to ferrite nomenclature
    intensity = FerriteIntensity(intensity_array, fgrid)

    #TODO: Check grid dimensions are compatible with start_img and spacing

    # instantiate generic grid
    FerriteImage{D,T,FSG}(intensity, num_pixels, start_img, spacing_img, fgrid)
end


"Creates a 2D FSG inside given a start point, `start_img`, spacing beteween pixels
`spacing_img`, length `length_img` and a number of pixels `num_pixels` in each direction."
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
    fgrid = FSG(start_grid_point, finish_grid_point, num_elements, element)

    return fgrid
end

"Creates a 3D FSG inside given a start point, `start_img`, spacing beteween pixels
`spacing_img`, length `length_img` and a number of pixels `num_pixels` in each direction."
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
    fgrid = FSG(start_grid_point, finish_grid_point, num_elements, element)

    return fgrid
end

"Returns an image intensity array of `img` in a ferrite grid.
A cooridnates transformation is done back to cartesian coordinates."
function intensity(img::AbstractImage{D,T,<:FSG}) where {D,T}

    intensity = img.intensity
    @assert typeof(intensity) <: FerriteIntensity
    intensity_vec = value(intensity)
    indexes = _indexes(intensity)

    intensity_array = Array{T,D}(undef, num_pixels(img))

    [intensity_array[indexes[i]] = int for (i, int) in enumerate(intensity_vec)]

    return intensity_array
end

"Interpolates the image `img` intensity in a vector of points `vec_points`."
function _interpolate(
    vec_points::Vector{NTuple{D,T}},
    img::AbstractImage{D,T,<:FSG},
) where {D,T}

    fgrid = grid(img).grid # ferrite grid

    intensity = value(img.intensity) # intensity with ferrite nomenclature

    # Scalar magnitude dofHandler
    dh = DofHandler(fgrid)
    push!(dh, :intensity, 1, Lagrange{D,RefCube,1}())
    close!(dh)

    # create a point evaluation handler
    eval_points = Vector{Vec{D,T}}()
    [push!(eval_points, Vec(p)) for p in vec_points]
    ph = PointEvalHandler(fgrid, eval_points)

    # evaluate magnitude at point p
    i_points = get_point_values(
        ph,
        dh,
        intensity,
        :intensity
    )

    return i_points

end

end# endmodule
