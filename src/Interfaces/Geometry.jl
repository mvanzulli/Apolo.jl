"""
Module defining geometric entities interface.
"""
module Geometry

import .Base: extrema, maximum, minimum

# Structured grid methods
export AbstractStructuredGrid
export ⊂, corners, coordinates, cartesian_index, dimension, extrema, element_type,
element_size, finish, grid, maximum, minimum, node_type, num_elements, num_nodes, start

# ==============
# Grids
# ==============
""" Abstract supertype for grids.

The following methods are provided by the interface:
- `⊂ (p,grid)`        -- returns true if a point is inside a grid
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
abstract type AbstractStructuredGrid{D,T} end


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
function grid(obj::Any)
    try
        obj.grid
    catch
        ERROR_GRID
    end
end

"Check if a point is inside a grid"
function ⊂(p::NTuple{D,T}, grid::AbstractStructuredGrid{D,T}) where {D,T}

    ex = extrema(grid)

    is_inside = true

    for axis in 1:D
        if is_inside == false
            return false
        else
            is_inside && ex[axis][D] ≤ p[D] ≤ ex[axis][D]
        end
    end
    return is_inside
end

"Interpolates a generic function inside a grid "
function _interpolate(
    p::NTuple{D},
    magnitude::M,
    grid::AbstractStructuredGrid{D}
) where {D,T,M}

    return ERROR_GRID
end

"Computes the dimension 2D Cartesian index given a 1D index in grid. "
function cartesian_index(onedim_index::Int, grid::AbstractStructuredGrid{2})

    num_nodes_grid = num_nodes(grid)

    (index_x, index_y) = _my_div_rem(onedim_index, num_nodes_grid[1], num_nodes_grid[1])

    return CartesianIndex(index_x, index_y)

end

"Computes the dimension 3D Cartesian index given a 1D index in grid. "
function cartesian_index(onedim_index::Int, grid::AbstractStructuredGrid{3})

    num_nodes_grid = num_nodes(grid)

    (index_z, rem_xy) = _my_div_rem(
        onedim_index,
        num_nodes_grid[1] * num_nodes_grid[2],
        num_nodes_grid[3]
        )

    (index_x, index_y) = _my_div_rem(rem_xy, num_nodes_grid[1], num_nodes_grid[1])

     return CartesianIndex(index_x, index_y, index_z)

end

"My own divisor reminder function"
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

end
