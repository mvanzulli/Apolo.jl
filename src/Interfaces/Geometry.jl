"""
Module defining geometric entities interface.
"""
module Geometry

import .Base: extrema, maximum, minimum

# Structured grid methods
export AbstractStructuredGrid
export ⊂, corners, coordinates, dimension, extrema, element_type, finish, grid, maximum, minimum,
    node_type, start

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
dimension(::AbstractStructuredGrid{D}) where D = D

"Returns the extrema of the grid."
extrema(::AbstractStructuredGrid) = error(ERROR_GRID)

"Returns the element type of the grid."
element_type(::AbstractStructuredGrid) = error(ERROR_GRID)

"Returns the grid final point."
finish(::AbstractStructuredGrid) = error(ERROR_GRID)

"Returns the maximum coordinate value at each direction."
maximum(::AbstractStructuredGrid) = error(ERROR_GRID)

"Returns the minimum coordinate value at each direction."
minimum(::AbstractStructuredGrid) = error(ERROR_GRID)

"Returns the node type"
node_type(::AbstractStructuredGrid{D,T}) where {D,T} = T

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
            is_inside  && ex[axis][D] ≤ p[D] ≤ ex[axis][D]
        end
    end
    return is_inside
end

end
