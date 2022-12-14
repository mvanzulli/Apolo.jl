
"""
Module Ferrite.jl grids implementations.
"""
module FerriteGrids

# TODO: Fix coordinates function overlead
import Apolo.Geometry: corners, coordinates, cartesian_index, dimension, extrema,
    element_type, element_size, finish, maximum, minimum, node_type, num_nodes, num_elements,
    start
import Apolo.Geometry: _extrapolate, _interpolate

using ..Geometry: AbstractStructuredGrid
using ..Geometry: _closest_point
using Ferrite: generate_grid, getcoordinates, getfaceset, getnodes, get_point_values, close!
using Ferrite: AbstractCell, CellIterator, DofHandler, FaceIndex, Grid, Lagrange,
    PointEvalHandler, Set, RefCube, Vec, Quadrilateral, Hexahedron

export FerriteStructuredGrid, border_points, set_coordinates

""" Creates a strcutred grid through ferrite grids.

## Fields
- `grid`     -- Ferrite.grid
- `vertices` -- Vertices vector
"""
struct FerriteStructuredGrid{D,E,T,V<:AbstractVector} <: AbstractStructuredGrid{D,T}
    grid::Grid{D,E,T}
    vertices::V
    num_elements::NTuple{D,<:Int}
end


# ==============
# Hard contracts
# ==============

grid(fgrid::FerriteStructuredGrid) = fgrid.grid

corners(fgrid::FerriteStructuredGrid) = fgrid.vertices

dimension(::FerriteStructuredGrid{D}) where {D} = D

function extrema(fgrid::FerriteStructuredGrid{D,E,T}) where {D,E,T}

    gird_corners = corners(fgrid)

    ext = Vector{NTuple{2,T}}()

    for axis in 1:D
        extrema_axis = extrema(getindex.(gird_corners, axis))
        push!(ext, extrema_axis)
    end

    return Tuple(ext)
end

element_type(::FerriteStructuredGrid{D,E}) where {D,E} = E

element_size(fgrid::FerriteStructuredGrid) = (finish(fgrid) .- start(fgrid)) ./ num_elements(fgrid)

node_type(::FerriteStructuredGrid{D,E,T}) where {D,E,T} = T

num_elements(fgrid::FerriteStructuredGrid) = fgrid.num_elements

num_nodes(fgrid::FerriteStructuredGrid{D}) where {D} = num_elements(fgrid) .+ Tuple(ones(Int, D))

maximum(fgrid::FerriteStructuredGrid) = maximum.(extrema(fgrid))

minimum(fgrid::FerriteStructuredGrid) = minimum.(extrema(fgrid))

start(fgrid::FerriteStructuredGrid) = getindex.(extrema(fgrid), 1)

finish(fgrid::FerriteStructuredGrid{D}) where {D} = getindex.(extrema(fgrid), 2)


" Creates a ferrite vector of corners in CCW "
function _corners(start_point::NTuple{2}, finish_point::NTuple{2})

    o???, o??? = start_point
    f???, f??? = finish_point


    grid_corners = [
        Vec{2}((o???, o???)),
        Vec{2}((f???, o???)),
        Vec{2}((f???, f???)),
        Vec{2}((o???, f???)),
    ]

    return grid_corners

end

function _corners(start_point::NTuple{3}, finish_point::NTuple{3})

    o???, o???, o??? = start_point
    f???, f???, f??? = finish_point

    grid_corners = [
        Vec{3}((o???, o???, o???)),
        Vec{3}((f???, o???, o???)),
        Vec{3}((f???, f???, o???)),
        Vec{3}((o???, f???, o???)),
        Vec{3}((o???, o???, f???)),
        Vec{3}((f???, o???, f???)),
        Vec{3}((f???, f???, f???)),
        Vec{3}((o???, f???, f???)),
    ]

    return grid_corners

end

"""
Creates a rectangular Ferrite grid with Quadrilateral.
"""
function _create_ferrite_rectangular_grid(
    start_point::NTuple{2},
    finish_point::NTuple{2},
    num_elements::NTuple{2},
    element::Type{ET}
) where {ET<:AbstractCell{2}}

    grid_corners = _corners(start_point, finish_point)

    return generate_grid(element, num_elements, grid_corners), grid_corners

end

"""
Creates a rectangular Ferrite grid with Hexahedron.
"""
function _create_ferrite_rectangular_grid(
    start_point::NTuple{3},
    finish_point::NTuple{3},
    num_elements::NTuple{3},
    element::Type{ET}
) where {ET<:Hexahedron}

    grid_corners = _corners(start_point, finish_point)

    # extract start and finish borders accordring to generate_grid ferrite function
    # the end -1 entry corresponds to the right ferrite nomenclature
    left = grid_corners[1]
    right = grid_corners[end-1]

    # return corners to add into FerriteGrid.vertices
    return generate_grid(element, num_elements, left, right), grid_corners

end

function FerriteStructuredGrid(
    start_point::NTuple{D},
    finish_point::NTuple{D},
    num_elements::NTuple{D},
    element::Type{ET}
) where {D,ET<:AbstractCell}

    grid_created, grid_corners = _create_ferrite_rectangular_grid(start_point, finish_point, num_elements, element)

    return FerriteStructuredGrid(grid_created, grid_corners, num_elements)

end

# ==============
# Features
# ==============
"DofHandler constructor using a FerriteStructuredGrid "
DofHandler(fgrid::FerriteStructuredGrid) = DofHandler(grid(fgrid))

"PointEvalHandler constructor using a FerriteStructuredGrid "
function PointEvalHandler(fgrid::FerriteStructuredGrid, eval_points)
    return PointEvalHandler(grid(fgrid), eval_points)
end

_borders(::FerriteStructuredGrid{2}) = ["top", "bottom", "left", "right"]
_borders(::FerriteStructuredGrid{3}) = ["top", "bottom", "left", "right", "back", "front"]

"Get ferrite facesets of the grid border"
function _boundary_indexes(fgrid::FerriteStructuredGrid)

    # extract Ferrite.Grid type
    ferrite_grid = grid(fgrid)

    ??faceset = Set{FaceIndex}()

    # get faces indexes
    for border in _borders(fgrid)

        faces_index = getfaceset(ferrite_grid, border)
        [push!(??faceset, face) for face in faces_index]

    end

    return ??faceset

end

" Gets a vector of Vec given a ferrite face set."
function set_coordinates(
    ??faceset::Set{FaceIndex},
    my_fgrid::FerriteStructuredGrid{D,E,T}
) where {D,E,T}

    # extract Ferrite.Grid type
    fgrid = grid(my_fgrid)

    # vector to fill
    ??faceset_coords = Vector{Vec{D,T}}()

    # iterate over each face and push the point if is not already included
    for face in ??faceset
        coords_face = getcoordinates(fgrid, face)
        for point in coords_face
            if (point ??? ??faceset_coords)
                push!(??faceset_coords, point)
            end
        end
    end

    return ??faceset_coords

end


"Returns the border points in vector of Vec given a ferrite grid. "
function border_points(
    fgrid::FerriteStructuredGrid{2},
    num_points_border::Int
)
    # check num_points_border value
    num_points_border < 2 && throw(ArgumentError("num_points_border must be greater than 2"))

    # build border ranges
    ((x?????????, x?????????), (y?????????, y?????????)) = extrema(fgrid)
    ??left_points = [Vec((x?????????, y)) for y in range(y?????????, y?????????, length=num_points_border)]
    ??top_points = [Vec((x, y?????????)) for x in range(x?????????, x?????????, length=num_points_border)]
    ??right_points = reverse([Vec((x?????????, y)) for y in range(y?????????, y?????????, length=num_points_border)])
    ??bottom_points = [Vec((x, y?????????)) for x in range(x?????????, x?????????, length=num_points_border)]


    # merge them into one single set
    merge_points(points, new_points) = [(p ??? points && push!(points, p)) for p in new_points]

    # initialize and merge
    ??grid_points = Vector{Vec{dimension(fgrid),node_type(fgrid)}}()
    merge_points(??grid_points, ??left_points)
    merge_points(??grid_points, ??top_points)
    merge_points(??grid_points, ??right_points)
    merge_points(??grid_points, ??bottom_points)

    return ??grid_points
end

" Build magnitude vector according to Ferrite's nomenclature (cell by cell CCW)"
function _convert_to_ferrite_nomenclature(
    mag::Array{T,D},
    fgrid::FerriteStructuredGrid{D}
) where {D,T<:Real}

    ferrite_magnitude = Vector{T}()

    # since cells share nodes (nodes that has been added are redundant)
    added_nodes_index = Vector{Int}()
    cartesian_index_nodes = Vector{CartesianIndex}()

    # iterate over each grid cell
    for cell in CellIterator(grid(fgrid))

        # get cell nodes
        nodes_cell = getnodes(cell)

        # add magnitude if the node has not been added yet
        for num_node_cell in nodes_cell

            if num_node_cell ??? added_nodes_index
                push!(added_nodes_index, num_node_cell)
                push!(
                    ferrite_magnitude,
                    mag[cartesian_index(num_node_cell, fgrid)]
                )
            end

        end

    end

    # cartesian indexes of each point
    cartesian_indexes = cartesian_index.(added_nodes_index, Ref(fgrid))

    return ferrite_magnitude, cartesian_indexes
end

#TODO: extend for vectorial magnitudes
"Interpolates a scalar magnitude with a ferrite grid"
function _interpolate(
    vec_points::Vector{NTuple{DG,T}},
    mag::Array{T,DM},
    mag_symbol::Symbol,
    fgrid::FerriteStructuredGrid{DG},
) where {DG,T,DM}

    [p ??? fgrid && throw(ArgumentError("p = $p ??? the grid domain")) for p in vec_points]
    DM > 3 && throw(ArgumentError("magnitude dimension cannot exceed 3"))

    # make the magnitude compatible with ferrite grids
    fmag, _ = _convert_to_ferrite_nomenclature(mag, fgrid)

    # Scalar magnitude dofHandler
    dh = DofHandler(fgrid)
    push!(dh, mag_symbol, 1, Lagrange{DG,RefCube,1}())
    close!(dh)

    # create a point evaluation handler
    eval_points = Vector{Vec{DG,T}}()
    [push!(eval_points, Vec(p)) for p in vec_points]
    ph = PointEvalHandler(fgrid, eval_points)

    # evaluate magnitude at point p
    m_points = get_point_values(
        ph,
        dh,
        fmag,
        mag_symbol
    )

    return m_points

end


"Interpolates a scalar magnitude with a ferrite grid"
function _extrapolate(
    vec_points::Vector{NTuple{DG,T}},
    mag::Array{T,DM},
    fgrid::FerriteStructuredGrid{DG}
) where {DG,T,DM}

    [p ??? fgrid &&
        throw(ArgumentError("p = $p ??? fgrid please use `_interpolate` method"))
     for p in vec_points]
    DM > 3 && throw(ArgumentError("magnitude dimension cannot exceed 3"))

    # make the magnitude compatible with ferrite grids
    fmag, _ = _convert_to_ferrite_nomenclature(mag, fgrid)

    dh = DofHandler(fgrid)
    push!(dh, :magnitude, 1, Lagrange{DG,RefCube,1}())
    close!(dh)

    # create a point evaluation handler
    eval_points = Vector{Vec{DG,T}}()

    [push!(eval_points, Vec(_closest_point(p, fgrid))) for p in vec_points]
    ph = PointEvalHandler(fgrid, eval_points)

    # evaluate magnitude at point p
    m_points = get_point_values(
        ph,
        dh,
        fmag,
        :magnitude
    )

    return m_points

end

end # endmodule
