"""
Module defining image properties and features.
"""
module Images

using Apolo.Geometry: AbstractStructuredGrid

import Base: length
import Apolo.Geometry: ⊂, ⊄, coordinates, cartesian_index, dimension, extrema, finish, grid, start
import Apolo.Geometry: _interpolate, _extrapolate

export AbstractImage
export finish_grid, intensity, intensity_type, num_pixels, start_grid, spacing


#################
# Abstract Image
#################
""" Abstract supertype for all images.

The following methods are provided by the interface:

- `coordinates(img)`       -- returns the image coordinates in three linear ranges.
- `dimension(img)`         -- returns the image dimension (1D, 2D or 3D).
- `extrema(img)`           -- returns a tuple of tuples containing max and min of each gird axis.
                            in the general case this should return
                            ((xₘᵢₙ, xₘₐₓ), (yₘᵢₙ, yₘₐₓ), (zₘᵢₙ, zₘₐₓ))
- `finish(img)`            -- returns the final point coordinates.
- `finish(img)`            -- returns the final point coordinates.
- `finish_grid(img)`       -- returns the final point of the image grid.
- `grid(img)`              -- returns the image grid.
- `intensity(img)`         -- returns the intensity array of the image.
- `intensity_type(img)`    -- returns the intensity data type.
- `length(img)`            -- returns the image dimensions in metric units.
- `num_pixels(img)`            -- returns the image resolution in pixels.
- `spacing(img)`           -- returns the space between pixels.
- `start(img)`             -- returns the start point coordinates.
- `start_grid(img)`        -- returns the start point of the image grid.

The following methods must be provided by the user:

_interpolate
_extrapolate

"""
abstract type AbstractImage{D,T,G} end

const ERROR_IMG = :("This method is not available for this image type. Please implement it")

" Gets the image coordinates"
function coordinates(img::AbstractImage)

    return LinRange.(start_grid(img), finish_grid(img), num_pixels(img))

end

cartesian_index(ondim_index::Int, img::AbstractImage) = cartesian_index(ondim_index, grid(img))

"Gets the image max and min coordinates"
function extrema(img::AbstractImage{D,T}) where {D,T}

    start_point = start(img)
    finish_point = finish(img)

    return [(start_point[axis], finish_point[axis]) for axis in 1:D]
end

" Gets the image dimension"
dimension(::AbstractImage{D,T}) where {D,T} = D

" Gets the image finish point"
finish(img::AbstractImage) = start(img) .+ length(img)
" Gets the image finish point"
function finish(
    start_img::NTuple{D,<:Real},
    num_pixels::NTuple{D,<:Integer},
    spacing_img::NTuple{D,<:Real},
) where {D}

    return start_img .+ num_pixels .* spacing_img
end
" Gets the final point"
function finish_grid(img::AbstractImage)
    try
        finish(grid(img))
    catch
        finish(img) .- spacing(img) ./ 2
    end
end

"Computes the grid final point."
function finish_grid(
    finish_img::NTuple{D,<:Real},
    spacing_img::NTuple{D,<:Real}
) where {D}
    return finish_img .- spacing_img ./ 2
end

" Gets the image grid"
function grid(img::AbstractImage)
    try
        img.grid
    catch
        error(ERROR_IMG)
    end
end

" Gets the image intensity array"
function intensity(img::AbstractImage)
    try
        img.intensity
    catch
        error(ERROR_IMG)
    end
end

" Gets the image intensity array data type"
intensity_type(::AbstractImage{D,T}) where {D,T} = T

" Gets the image metric dimensions"
function length(img::AbstractImage)
    try
        num_pixels(img) .* spacing(img)
    catch
        error(ERROR_IMG)
    end
end

function length(
    start_img::NTuple{D,<:Real},
    finish_img::NTuple{D,<:Real}
) where {D}
    return finish_img .- start_img
end

" Gets the image dimensions in voxels"
function num_pixels(img::AbstractImage)
    try
        img.num_pixels
    catch
        error(ERROR_IMG)
    end
end

"Computes the grid number of elements."
function _grid_num_elements(
    num_pixels::NTuple{D,<:Integer}
) where {D}
    return num_pixels .- Tuple(ones(Int, D))
end

" Gets the image spacing between pixels"
function spacing(img::AbstractImage)
    try
        img.spacing
    catch
        error(ERROR_IMG)
    end
end

" Gets the image start point"
function start(img::AbstractImage)
    try
        img.start
    catch
        error(ERROR_IMG)
    end
end

" Gets the image grid origin"
function start_grid(img::AbstractImage)
    try
        start(grid(img))
    catch
        start(img) .+ spacing(img) ./ 2
    end
end

"Computes the grid start point."
function start_grid(
    start_img::NTuple{D,<:Real},
    spacing_img::NTuple{D,<:Real}
) where {D}
    return start_img .+ spacing_img ./ 2
end

" Gets the total number of pixels"
_total_num_pixels(img::AbstractImage) = prod(num_pixels(img))

"Checks if a index is inbounds the image dimensions"
_index_is_inbounds(onedim_index::Int, img::AbstractImage) = 1 ≤ onedim_index ≤ _total_num_pixels(img)

"Checks if a index is inbounds the image dimensions"
function _index_is_inbounds(ndim_index::CartesianIndex{D}, img::AbstractImage{D,T}) where {D,T}

    num_pixels = num_pixels(img)

    # Check all the indexes are inbounds
    index_is_in = true
    for direction in 1:D
        index_is_in = index_is_in && 1 ≤ ndim_index[direction] ≤ num_pixels[direction]
    end

    return index_is_in

end

"Checks if a point is inbounds the image "
function ⊂(p::NTuple{D}, img::AbstractImage{D}) where {D}

    ex = extrema(img)
    for axis in 1:D

        is_in = ex[axis][1] ≤ p[axis] ≤ ex[axis][end]
        !is_in && return false

    end
    return true

end

"Checks if a point is outside a grid"
⊄(p::NTuple, img::AbstractImage) = !(⊂(p, img))

(img::AbstractImage{2} where {T})(x, y) = _eval_intensity((x, y), img)
(img::AbstractImage{3} where {T})(x, y, z) = print("IMPLEMENT me")

"Internal function to evaluate the image intensity at a generic point"
function _eval_intensity(
    p::NTuple{D,T},
    img::AbstractImage{D}
) where {D,T}

    # Check p is inisde the image
    p ⊂ img ? nothing : throw(BoundsError("p = $p is not inside the img frame"))
    # Check if p is inside the grid image
    if p ⊂ grid(img)
        intensity_p = _interpolate([p], img)
    else
        intensity_p = _extrapolate([p], img)
    end
    return getindex(intensity_p)
end

"Interpolates the image intensity in a vector of points "
function _interpolate(
    vec_points::Vector{NTuple{D,T}},
    img::AbstractImage{D}
) where {D,T}

    grid_img = grid(img)

    intensity_img = intensity(img)

    return _interpolate(vec_points, intensity_img, grid_img)
end

"Extrapolates the image intensity inside a vector of points "
function _extrapolate(
    vec_points::Vector{NTuple{D,T}},
    img::AbstractImage{D}
) where {D,T}

    grid_img = grid(img)

    intensity_img = intensity(img)

    return _extrapolate(vec_points, intensity_img, grid_img)
end

####################
# Abstract Intensity
####################
""" Abstract supertype for image intensities.

The following methods are provided by the interface:

- `dimension(intensity)`   -- returns the intensity array dimension.
- `value(intensity)`       -- returns the image intensity array.
- `typeof(intensity)`      -- returns the intensity type.
- `maximum(intensity)`     -- returns the intensity maximum value.
- `minimum(intensity)`     -- returns the intensity maximum value.
"""
abstract type AbstractIntensity{D,T} end

" Gets the intensity array value"
function value(int::AbstractIntensity)
    try
        int.value
    catch
        error(ERROR_IMG)
    end
end

intensity_type(::AbstractIntensity{D,T}) where {D,T} = T

dimension(::AbstractIntensity{D,T}) where {D,T} = D

# ====================
# Images implementations
# ====================
include("../Images/ferrite_img.jl")


# ####################
# # Generic Image
# ####################

# """ Generic image struct.

# ### Fields:

# - `intensity`  -- Intensity array
# - `num_pixels` -- Number of pixels in each direction
# - `spacing`    -- Space in each direction between each pixel
# - `start`      -- Start coordinates (considering the start located at [1,1,1] )
# - `grid`       -- Ferrite grid inside the image
# """
# struct GenericImage{T,D} <: AbstractImage{T,D}
#     intensity::Array{T,D}
#     num_pixels::NTuple{D,<:Integer}
#     spacing::NTuple{D,<:Real}
#     start::NTuple{D,<:Real}
#     grid::Grid

#     # Cash the Ferrite.Grid into grid field of MedicalImage
#     function GenericImage(
#         intensity::Array{T,D},
#         num_pixels::NTuple{D,<:Integer},
#         spacing::NTuple{D,<:Real},
#         start::NTuple{D,<:Real},
#     ) where {T,D}

#         # compute end point
#         finish = start .+ num_pixels .* spacing
#         length_img = finish .- start
#         # create grid
#         grid = create_img_grid(start, spacing, length_img, num_pixels)

#         # instantiate generic grid
#         new{T,D}(intensity, num_pixels, spacing, start, grid)
#     end
# end




end
