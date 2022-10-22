"""
Module defining geometric entities interface.
"""
module Geometry

import .Base: extrema, maximum, minimum

# Structured grid methods
export corners, dimension, extrema, finish, maximum, minimum, node_type, start

# ==============
# Grids
# ==============
""" Abstract supertype for grids.

The following methods are provided by the interface:

- `start(grid)`       -- returns the grid start point.
- `finish(grid)`      -- returns the grid final point.
- `maximum(grid)`     -- returns a tuple with the maximum coordinate value of the
                        grid towards an axis.
- `minimum(grid)`     -- returns a tuple with the minimum coordinate value of the
                        grid towards an axis.
- `extrema(grid)`     -- returns a tuple of tuples containing max and min of each gird axis.
- `corners(grid)`     -- returns the grid coordinates.
- `boundary(grid)`    -- returns the grid boundaries.
- `coordinates(grid)` -- returns the grid coordinates.
- `dimension(grid)`   -- returns the grid dimension.
- `node_type(grid)`   -- returns the type of each grid node.
"""
abstract type AbstractStructuredGrid{D,T} end


const ERROR_message = :("This grid type has not this method implemented ")

"Returns the grid start point."
start(::AbstractStructuredGrid) = error(ERROR_message)

"Returns the grid final point."
finish(::AbstractStructuredGrid) = error(ERROR_message)

"Returns the extrema of the grid."
extrema(::AbstractStructuredGrid) = error(ERROR_message)

"Returns the maximum coordinate value at each direction."
maximum(::AbstractStructuredGrid) = error(ERROR_message)

"Returns the minimum coordinate value at each direction."
minimum(::AbstractStructuredGrid) = error(ERROR_message)

"Returns the grid corners"
corners(::AbstractStructuredGrid) = error(ERROR_message)

"Returns the grid boundary"
boundary(::AbstractStructuredGrid) = error(ERROR_message)

"Returns the grid coordinates"
coordinates(::AbstractStructuredGrid) = error(ERROR_message)

"Returns the grid dimension"
dimension(::AbstractStructuredGrid{D}) where D = D

"Returns the node type"
node_type(::AbstractStructuredGrid{D,T}) where {D,T} = T


end
