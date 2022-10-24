# Prototyping
using Test: @test
using Apolo
# using Apolo: _corners, _boundary_indexes
using Ferrite: Quadrilateral, Grid, Set, FaceIndex, Vec

const INTERVAL_START = LinRange(-10.,10.,20)

include("ferrite_grids.jl")

# Define a random grid
start_point = (rand(INTERVAL_START), rand(INTERVAL_START))
finish_point = start_point .+ (rand(.1:10.), rand(.1:10.))
num_elements_grid = (rand(1:5), rand(1:5))
elemtype = Quadrilateral
fgrid = FerriteStructuredGrid(start_point, finish_point, num_elements_grid, elemtype)

# test getter functions
@test corners(fgrid) == _corners(start_point, finish_point)
@test extrema(fgrid) == ((start_point[1], finish_point[1]), (start_point[end], finish_point[end]))
@test dimension(fgrid) == length(start_point)
@test element_type(fgrid) == elemtype
@test node_type(fgrid) == eltype(finish_point)
@test start(fgrid) == start_point
@test finish(fgrid) == finish_point


# Iterate over all fields of the object and check if are equal
grid_to_test = grid(fgrid)
fgrid_fields = fieldnames(typeof(grid_to_test))

grid_test_bench,_ = _create_ferrite_rectangular_grid(
    start_point,
    finish_point,
    num_elements_grid,
    elemtype
)

grids_are_equal = true
for field in fgrid_fields
    grids_are_equal &&
    getproperty(grid_test_bench, field) == getproperty(grid_to_test, field)
end

@test start_point ⊂ fgrid && finish_point ⊂ fgrid &&
 (finish_point .+ start_point) ./ 2 ⊂ ferrite_grid


# test more specific ferrite functions
Ωsets = _boundary_indexes(fgrid)

# all sets must be unique
@test length(unique(Ωsets)) == length(Ωsets)

# test start point is in the first border cell
start_cell = Set{FaceIndex}()
push!(start_cell, FaceIndex((1, 1)))
start_cell_coords = coordinates(start_cell, fgrid)
@test Vec{dimension(fgrid),node_type(fgrid)}(start_point) ∈ start_cell_coords

# test corners are in border points
points_per_segment = rand(2:5)
Ωfgrid_points = border_points(fgrid, points_per_segment)
bool_corners_in_borders = true
for point in corners(fgrid)
    bool_corners_in_borders && point ∈ Ωfgrid_points
end
@test bool_corners_in_borders
