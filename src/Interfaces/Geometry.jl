"""
Module defining geometric entities interface.
"""
module Geometry

import .Base: extrema, maximum, minimum

# Structured grid methods
export AbstractStructuredGrid
export corners, coordinates, cartesian_index, dimension, extrema, element_type,
    element_size, finish, grid, maximum, minimum, node_type, num_elements, num_nodes, start

# ==============
# Abstract Grid
# ==============
""" Abstract supertype for structured grids.

The following methods are provided by the interface:
- `∈ (p,grid)`        -- returns true if a point is inside a grid.
- `∉ (p,grid)`        -- returns true if a point is outside a grid.
- `boundary(grid)`    -- returns the grid boundaries.
- `corners(grid)`     -- returns the grid corners.
- `coordinates(grid)` -- returns the grid coordinates.
- `dimension(grid)`   -- returns the grid dimension.
- `element_type(grid)`-- returns the element type of the grid.
- `element_size(grid)`-- returns the element size of the grid.
- `extrema(grid)`     -- returns a tuple of tuples containing max and min of each gird axis.
                         in the general case this should return
                         ((xₘᵢₙ, xₘₐₓ), (yₘᵢₙ, yₘₐₓ), (zₘᵢₙ, zₘₐₓ))
- `finish(grid)`      -- returns the grid final point.
- `grid(object)`      -- returns the grid of an object.
- `maximum(grid)`     -- returns a tuple with the maximum coordinate value of the grid
                         towards an axis.
- `minimum(grid)`     -- returns a tuple with the minimum coordinate value of the grid
                         towards an axis.
- `node_type(grid)`   -- returns the type of each grid node.
- `num_nodes(grid)`   -- returns the number of grid nodes.
- `start(grid)`       -- returns the grid start point.
"""
abstract type AbstractGrid{D,T} end

abstract type AbstractStructuredGrid{D,T} <: AbstractGrid{D,T} end


const ERROR_GRID = :("This method is not available for this grid type. Please implement it")

"Returns the grid boundary"
boundary(::AbstractStructuredGrid) = error(ERROR_GRID)

"Returns the grid corners"
corners(::AbstractStructuredGrid) = error(ERROR_GRID)

"Returns the grid coordinates"
coordinates(::AbstractStructuredGrid) = error(ERROR_GRID)

"Returns the grid dimension"
dimension(::AbstractStructuredGrid{D}) where {D} = D

"Returns the extrema of the grid."
extrema(::AbstractStructuredGrid) = error(ERROR_GRID)

"Returns the element type of the grid."
element_type(::AbstractStructuredGrid) = error(ERROR_GRID)

"Returns the element element size in x,y of the grid."
element_size(::AbstractStructuredGrid) = error(ERROR_GRID)

"Returns the grid final point."
finish(::AbstractStructuredGrid) = error(ERROR_GRID)

"Returns the maximum coordinate value at each direction."
maximum(::AbstractStructuredGrid) = error(ERROR_GRID)

"Returns the minimum coordinate value at each direction."
minimum(::AbstractStructuredGrid) = error(ERROR_GRID)

"Returns the node type"
node_type(::AbstractStructuredGrid{D,T}) where {D,T} = T

"Returns a tuple with the number of nodes in each direction"
num_nodes(::AbstractStructuredGrid) = error(ERROR_GRID)

"Returns a tuple with the number of elements in each direction"
num_elements(::AbstractStructuredGrid) = error(ERROR_GRID)

"Returns the grid start point."
start(::AbstractStructuredGrid) = error(ERROR_GRID)

"Returns the grid of an object"
grid(obj::Any) = obj.grid

"Checks if a point is inside a grid"
function Base.:∈(p::NTuple{D,T}, grid::AbstractStructuredGrid{D,T}) where {D,T}

    ex = extrema(grid)

    for axis in 1:D

        is_in = ex[axis][1] ≤ p[axis] ≤ ex[axis][end]
        !is_in && return false

    end
    return true
end

"Interpolates a generic function inside a grid "
function _interpolate(
    vec_points::Vector{NTuple{D,T}},
    magnitude::M,
    mag_symbol::Symbol,
    grid::AbstractStructuredGrid{D},
) where {D,T,M}

    return ERROR_GRID
end

"Computes the dimension 2D Cartesian index given a 1D index in grid. "
function cartesian_index(onedim_index::Int, grid::AbstractStructuredGrid{2})

    onedim_index ≤ 0 && throw(ArgumentError("onedim_index = $onedim_index ≤ 0"))

    num_nodes_grid = num_nodes(grid)

    (index_y, rem_x) = divrem(onedim_index, num_nodes_grid[1])
    if rem_x == 0
        index_x = num_nodes_grid[1]
        return CartesianIndex(index_x, index_y)
    end
    return CartesianIndex(rem_x, index_y + 1)
end

"Computes the dimension 3D Cartesian index given a 1D index in grid. "
function cartesian_index(onedim_index::Int, grid::AbstractStructuredGrid{3})

    onedim_index ≤ 0 && throw(ArgumentError("onedim_index = $onedim_index ≤ 0"))

    num_nodes_grid = num_nodes(grid)

    (index_z, rem_xy) = divrem(onedim_index, num_nodes_grid[1] * num_nodes_grid[2])

    if rem_xy == 0
        index_x, index_y = num_nodes_grid[1:2]
        return CartesianIndex(index_x, index_y, index_z)
    else
        (index_y, rem_x) = divrem(rem_xy, num_nodes_grid[1])
        if rem_x == 0
            index_x = num_nodes_grid[1]
            return CartesianIndex(index_x, index_y, index_z + 1)
        end
        return CartesianIndex(rem_x, index_y + 1, index_z + 1)
    end

end

"Computes which grid border should be used to extrapolate for a 2D grid"
function _which_border(
    p::NTuple{2},
    start_point::NTuple{2},
    finish_point::NTuple{2},
)

    # Bottom borders
    #(1,1)
    p[1] ≤ start_point[1] && p[2] ≤ start_point[2] && return :left_bottom
    #(end,1)
    p[1] ≥ finish_point[1] && p[2] ≤ start_point[2] && return :right_bottom

    # Top borders
    #(1,end)
    p[1] ≤ start_point[1] && p[2] ≥ finish_point[2] && return :left_top
    #(end,end)
    p[1] ≥ finish_point[1] && p[2] ≥ finish_point[2] && return :right_top

    # Simple cases
    p[2] < start_point[2] && return :bottom
    p[1] < start_point[1] && return :left
    p[2] > finish_point[2] && return :top
    p[1] > finish_point[1] && return :right

end

"Computes which grid border should be used to extrapolate for 3D grids"
function _which_border(
    p::NTuple{3},
    start_point::NTuple{3},
    finish_point::NTuple{3},
)

    # 2D classification of the point
    twodim_border = _which_border(
        (p[1], p[2]),
        (start_point[1], start_point[2]),
        (finish_point[1], finish_point[2]),
    )

    if p[3] ≤ start_point[3] # Front cases

        # Bottom borders
        twodim_border == :left_bottom && return :left_bottom_front
        twodim_border == :right_bottom && return :right_bottom_front

        # Top borders
        twodim_border == :left_top && return :left_top_front
        twodim_border == :right_top && return :right_top_front

        # Simple cases
        twodim_border == :bottom && return :bottom_front
        twodim_border == :left && return :left_front
        twodim_border == :top && return :top_front
        twodim_border == :right && return :right_front

        # If is in the previous regions then just front
        return :front

    elseif p[3] ≥ finish_point[3] # Back cases

        # Bottom borders
        twodim_border == :left_bottom && return :left_bottom_back
        twodim_border == :right_bottom && return :right_bottom_back

        # Top borders
        twodim_border == :left_top && return :left_top_back
        twodim_border == :right_top && return :right_top_back

        # Simple cases
        twodim_border == :bottom && return :bottom_back
        twodim_border == :left && return :left_back
        twodim_border == :top && return :top_back
        twodim_border == :right && return :right_back

        # If is in the previous regions then just back
        return :back

    elseif start_point[3] < p[3] < finish_point[3] # 2D case with zₘᵢₙ <  z < zₘₐₓ

        return twodim_border

    end

end


"Computes the border closest point to extrapolate considering 2D grids"
function _closest_point(
    p::NTuple{2},
    fgrid::AbstractStructuredGrid{2}
)

    p ∈ fgrid && throw(ArgumentError("p = $p ∈ fgrid please use `_interpolate` method"))

    start_point = start(fgrid)
    finish_point = finish(fgrid)

    border = _which_border(p, start_point, finish_point)

    border == :left_bottom && return start_point
    border == :left_top && return (start_point[1], finish_point[2])
    border == :right_top && return finish_point
    border == :right_bottom && return (finish_point[1], start_point[2])
    border == :bottom && return (p[1], start_point[2])
    border == :top && return (p[1], finish_point[2])
    border == :right && return (finish_point[1], p[2])
    border == :left && return (start_point[1], p[2])

end

"Computes the border closest point to extrapolate considering 3D grids"
function _closest_point(
    p::NTuple{3},
    fgrid::AbstractStructuredGrid{3},
)

    p ∈ fgrid && throw(ArgumentError("p = $p ∈ fgrid please use `_interpolate` method"))

    start_point = start(fgrid)
    finish_point = finish(fgrid)

    border = _which_border(p, start_point, finish_point)

    # Front borders
    border == :left_bottom_front && return start_point
    border == :right_bottom_front && return (finish_point[1], start_point[2], start_point[3])
    border == :left_top_front && return (start_point[1], finish_point[2], start_point[3])
    border == :right_top_front && return (finish_point[1], finish_point[2], start_point[3])
    border == :bottom_front && return (p[1], start_point[2], start_point[3])
    border == :top_front && return (p[1], finish_point[2], start_point[3])
    border == :left_front && return (start_point[1], p[2], start_point[3])
    border == :right_front && return (finish_point[1], p[2], start_point[3])
    border == :front && return (p[1], p[2], start_point[3])

    # Back borders
    border == :left_bottom_back && return (start_point[1], start_point[2], finish_point[3])
    border == :right_bottom_back && return (finish_point[1], start_point[2], finish_point[3])
    border == :left_top_back && return (start_point[1], finish_point[2], finish_point[3])
    border == :right_top_back && return finish_point
    border == :bottom_back && return (p[1], start_point[2], finish_point[3])
    border == :top_back && return (p[1], finish_point[2], finish_point[3])
    border == :left_back && return (start_point[1], p[2], finish_point[3])
    border == :right_back && return (finish_point[1], p[2], finish_point[3])
    border == :back && return (p[1], p[2], finish_point[3])

    # zₘᵢₙ <  z < zₘₐₓ borders
    border == :bottom && return (p[1], start_point[2], p[3])
    border == :top && return (p[1], finish_point[2], p[3])
    border == :right && return (finish_point[1], p[2], p[3])
    border == :left && return (start_point[1], p[2], p[3])


    border == :left_bottom && return (start_point[1], start_point[2], p[3])
    border == :left_top && return (start_point[1], finish_point[2], p[3])
    border == :right_top && return (finish_point[1], finish_point[2], p[3])
    border == :right_bottom && return (finish_point[1], start_point[2], p[3])

end


# ====================
# Grids implementations
# ====================

include("../Geometry/ferrite_grids.jl")

end
