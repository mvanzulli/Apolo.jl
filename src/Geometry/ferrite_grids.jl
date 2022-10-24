
##########################################################
# Main types and functions to handle with ferrite grids #
##########################################################

# TODO: Fix coordinates function overlead
import Apolo.Geometry: corners, coordinates, dimension, extrema, element_type, finish,
    maximum, minimum, node_type, start

using Apolo.Geometry: AbstractStructuredGrid
using Ferrite: Grid, generate_grid, getcoordinates, getfaceset
using Ferrite: AbstractCell, FaceIndex, Set, Vec

export FerriteStructuredGrid, border_points, set_coordinates

""" Creates a strcutred grid through ferrite grids.

## Fields
- `grid`     -- Ferrite.grid
- `vertices` -- Vertices vector
"""
struct FerriteStructuredGrid{D,E,T} <: AbstractStructuredGrid{D,T}
    grid::Grid{D,E,T}
    vertices::AbstractVector
end

grid(fgrid::FerriteStructuredGrid) = fgrid.grid

corners(fgrid::FerriteStructuredGrid) = fgrid.vertices

dimension(::FerriteStructuredGrid{D}) where {D} = D

function extrema(fgrid::FerriteStructuredGrid{D,E,T}) where {D,E,T}

    gird_corners = corners(fgrid)

    ext = Vector{NTuple{D,T}}()

    for axis in 1:D
        extrema_axis = extrema(getindex.(gird_corners, axis))
        push!(ext, extrema_axis)
    end

    return Tuple(ext)
end

element_type(::FerriteStructuredGrid{D,E}) where {D,E} = E

node_type(::FerriteStructuredGrid{D,E,T}) where {D,E,T} = T

maximum(fgrid::FerriteStructuredGrid) = maximum.(extrema(fgrid))

minimum(fgrid::FerriteStructuredGrid) = minimum.(extrema(fgrid))

start(fgrid::FerriteStructuredGrid) = getindex.(extrema(fgrid), 1)

finish(fgrid::FerriteStructuredGrid{D}) where {D} = getindex.(extrema(fgrid), D)


" Creates a ferrite vector of corners in CCW "
function _corners(start_point::NTuple{2}, finish_point::NTuple{2})

    oᵢ, oⱼ = start_point
    fᵢ, fⱼ = finish_point


    grid_corners = [
        Vec{2}((oᵢ, oⱼ)),
        Vec{2}((fᵢ, oⱼ)),
        Vec{2}((fᵢ, fⱼ)),
        Vec{2}((oᵢ, fⱼ)),
    ]

    return grid_corners

end

function _corners(start_point::NTuple{3}, finish_point::NTuple{3})

    oᵢ, oⱼ, oₖ = start_point
    fᵢ, fⱼ, fₖ = finish_point

    grid_corners = [
        Vec{3}((oᵢ, oⱼ, oₖ)),
        Vec{3}((fᵢ, oⱼ, oₖ)),
        Vec{3}((fᵢ, fⱼ, oₖ)),
        Vec{3}((oᵢ, fⱼ, oₖ)),
        Vec{3}((oᵢ, oⱼ, fₖ)),
        Vec{3}((fᵢ, oⱼ, fₖ)),
        Vec{3}((fᵢ, fⱼ, fₖ)),
        Vec{3}((oᵢ, fⱼ, fₖ)),
    ]

    return grid_corners

end

"""
Creates a rectangular Ferrite grid.
"""
function _create_ferrite_rectangular_grid(
    start_point::NTuple{D},
    finish_point::NTuple{D},
    num_elements::NTuple{D},
    element::Type{ET}
) where {D,ET<:AbstractCell}

    grid_corners = _corners(start_point, finish_point)

    return generate_grid(element, num_elements, grid_corners), grid_corners
end

function FerriteStructuredGrid(
    num_elements::NTuple{D},
    start_point::NTuple{D},
    finish_point::NTuple{D},
    element::Type{ET}
) where {D,ET<:AbstractCell}

    grid_created, grid_corners = _create_ferrite_rectangular_grid(num_elements, start_point, finish_point, element)

    return FerriteStructuredGrid(grid_created, grid_corners)

end


"Get ferrite facesets of the grid border"
function _boundary_indexes(my_fgrid::FerriteStructuredGrid{2})

    # extract Ferrite.Grid type
    fgrid = grid(my_fgrid)

    Ωfaceset = Set{FaceIndex}()

    # get faces indexes
    for border in ["top", "bottom", "left", "right"]

        faces_index = getfaceset(fgrid, border)
        [push!(Ωfaceset, face) for face in faces_index]

    end

    return Ωfaceset

end

" Gets a vector of Vec given a ferrite face set."
function set_coordinates(
    Ωfaceset::Set{FaceIndex},
    my_fgrid::FerriteStructuredGrid{D,E,T}
) where {D,E,T}

    # extract Ferrite.Grid type
    fgrid = grid(my_fgrid)

    # vector to fill
    Ωfaceset_coords = Vector{Vec{D,T}}()

    # iterate over each face and push the point if is not already included
    for face in Ωfaceset
        coords_face = getcoordinates(fgrid, face)
        for point in coords_face
            if (point ∉ Ωfaceset_coords)
                push!(Ωfaceset_coords, point)
            end
        end
    end

    return Ωfaceset_coords

end


"Returns the border points in vector of Vec given a ferrite grid. "
function border_points(
    fgrid::FerriteStructuredGrid{2},
    num_points_border::Int
)
    # check num_points_border value
    num_points_border < 2 && throw(ArgumentError("num_points_border must be greater than 2"))

    # build border ranges
    ((xₘᵢₙ, xₘₐₓ), (yₘᵢₙ, yₘₐₓ)) = extrema(fgrid)
    Ωleft_points = [Vec((xₘᵢₙ, y)) for y in range(yₘᵢₙ, yₘₐₓ, length=num_points_border)]
    Ωtop_points = [Vec((x, yₘₐₓ)) for x in range(xₘᵢₙ, xₘₐₓ, length=num_points_border)]
    Ωright_points = reverse([Vec((xₘₐₓ, y)) for y in range(yₘᵢₙ, yₘₐₓ, length=num_points_border)])
    Ωbottom_points = [Vec((x, yₘᵢₙ)) for x in range(xₘᵢₙ, xₘₐₓ, length=num_points_border)]


    # merge them into one single set
    merge_points(points, new_points) = [(p ∉ points && push!(points, p)) for p in new_points]

    # initialize and merge
    Ωgrid_points = Vector{Vec{dimension(fgrid),node_type(fgrid)}}()
    merge_points(Ωgrid_points, Ωleft_points)
    merge_points(Ωgrid_points, Ωtop_points)
    merge_points(Ωgrid_points, Ωright_points)
    merge_points(Ωgrid_points, Ωbottom_points)

    return Ωgrid_points
end
