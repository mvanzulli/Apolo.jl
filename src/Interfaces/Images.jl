"""
Module defining image properties and features.
"""
module Images

using ..Geometry: AbstractStructuredGrid
using Reexport: @reexport

@reexport import ..Geometry:  ⊂, ⊄, coordinates, cartesian_index, dimension, extrema, finish,
grid, start
import ..Geometry: _interpolate, _extrapolate

export AbstractImage
export finish_grid, intensity, intensity_type, num_pixels, start_grid, spacing, path, value


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
"""
abstract type AbstractImage{D,T,G} end

const ERROR_IMG = :("This method is not available for this image type. If corresponds please implement it")

" Gets the image coordinates"
function coordinates(img::AbstractImage)

    return LinRange.(start_grid(img), finish_grid(img), num_pixels(img))

end

cartesian_index(ondim_index::Int, img::AbstractImage) = cartesian_index(ondim_index, grid(img))

"Gets the image max and min coordinates"
function extrema(img::AbstractImage{D}) where {D}

    start_point = start(img)
    finish_point = finish(img)

    return [(start_point[axis], finish_point[axis]) for axis in 1:D]
end

" Gets the image dimension"
dimension(::AbstractImage{D,T}) where {D,T} = D

" Gets the image finish point"
function finish(img::AbstractImage)
    try
        return img.finish
    catch
        return  start(img) .+ length(img)
    end
end

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
function Base.length(img::AbstractImage)
        try
            num_pixels(img) .* spacing(img)
        catch
            error(ERROR_IMG)
        end
end

function Base.length(
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

"Gets the image path if corresponds "
function path(img::AbstractImage)
    try
        img.path
    catch
        error(ERROR_IMG)
    end
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

"Evaluation functor for a 2D AbstractImage"
function (img::AbstractImage{2})(x::T, y::T; offset::NTuple{2,T} = Tuple(zeros(T,2))) where {T}
    _eval_intensity((x, y), img, offset = offset)
end

"Evaluation functor for a 3D AbstractImage"
function (img::AbstractImage{3})(x::T, y::T, z::T; offset::NTuple{3,T} = Tuple(zeros(T,3))) where {T}
    _eval_intensity((x, y, z), img, offset = offset)
end

"Evaluation functor for a vector of points"
function (img::AbstractImage{D})(
    vec_points::Vector{NTuple{D,T}};
    offset::NTuple{D,T} = Tuple(zeros(T,D))
) where {D,T}
    return _eval_intensity(vec_points, img, offset = offset)
end

"Internal function to evaluate the image intensity at a generic point"
function _eval_intensity(
    p::NTuple{D,T},
    img::AbstractImage{D};
    offset::NTuple{D,T} = Tuple(zeros(T,D))
) where {D,T}

    p = p .+ offset

    # Check p is inisde the image
    p ⊂ img ? nothing : throw(ArgumentError("p + offset = $p is not inside the img frame"))

    # Check if p is inside the grid image
    if p ⊂ grid(img)
        intensity_p = _interpolate([p], img)
    else
        intensity_p = _extrapolate([p], img)
    end
    return getindex(intensity_p)
end

"Internal function to evaluate the image intensity at a vector of points"
function _eval_intensity(
    vec_points::Vector{NTuple{D,T}},
    img::AbstractImage{D};
    offset::NTuple{D,T} = Tuple(zeros(T,D))
) where {D,T}

    intensity_vec = Vector{T}(undef, length(vec_points))

    # Check if p is inside the grid image
    for (num_point, p) in enumerate(vec_points)
        p = p .+ offset
        if p ⊂ img
            if p ⊂ grid(img)
                intensity_p = _interpolate([p], img)
            else
                intensity_p = _extrapolate([p], img)
            end
            intensity_vec[num_point] = getindex(intensity_p)
        else
            return missing
            throw(Warning("p + offset = $p is not inside the img frame"))
        end
    end

    return intensity_vec

end

"Interpolates the image intensity in a vector of points "
function _interpolate(vec_points::Vector{NTuple{D,T}},img::AbstractImage{D}) where {D,T}

    grid_img = grid(img)

    intensity_img = intensity(img)

    return _interpolate(vec_points, intensity_img, :intensity, grid_img)
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

- `value(intensity)`       -- returns the image intensity array.
- `intensity_type(intensity)`      -- returns the intensity type.
- `maximum(intensity)`     -- returns the intensity maximum value.
- `minimum(intensity)`     -- returns the intensity minimum value.
- `extrema(intensity)`     -- returns max and min intensity values.
"""
abstract type AbstractIntensity{T} end

" Gets the intensity array value"
function value(int::AbstractIntensity)
    try
        int.value
    catch
        error(ERROR_IMG)
    end
end

"Computes the maximum intensity"
maximum(I::AbstractIntensity) = maximum(value(I))

"Computes the mimunum intensity"
minimum(I::AbstractIntensity) = minimum(value(I))

"Computes the intensity type"
intensity_type(::AbstractIntensity{T}) where {T} = T

struct ScalarIntensity{T} <: AbstractIntensity{T}
    value::Array{T}
end


# ====================
# Images implementations
# ====================
include("../Images/ferrite_img.jl")
include("../Images/analytic_img.jl")
include("../Images/vtk_img.jl")
include("../Images/medical_img.jl")




end
