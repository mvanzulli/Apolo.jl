##############
# Grid tests #
##############

using Test: @test, @testset
using Apolo.Geometry
using Apolo: FerriteStructuredGrid, border_points, coordinates,
_convert_to_ferrite_nomenclature, _create_ferrite_rectangular_grid,
_corners, _boundary_indexes, _interpolate

using Ferrite: Quadrilateral, Grid, Set, FaceIndex, Vec

const INTERVAL_START = LinRange(-10.0, 10.0, 20)

@testset "Ferrite grids unitary tests" begin

    # Define a random grid
    start_point = (rand(INTERVAL_START), rand(INTERVAL_START))
    start_point = (1., 1.)
    finish_point = start_point .+ (rand(0.1:10.0), rand(0.1:10.0))
    finish_point = (2., 2.)
    num_elements_grid = (2, 2)
    elemtype = Quadrilateral
    fgrid = FerriteStructuredGrid(start_point, finish_point, num_elements_grid, elemtype)

    # test getter functions
    @test corners(fgrid) == _corners(start_point, finish_point)
    @test extrema(fgrid) == ((start_point[1], finish_point[1]), (start_point[end], finish_point[end]))
    @test dimension(fgrid) == length(start_point)
    @test element_size(fgrid) == (finish_point .- start_point) ./ num_elements_grid
    @test element_type(fgrid) == elemtype
    @test node_type(fgrid) == eltype(finish_point)
    @test num_elements(fgrid) == num_elements_grid
    @test num_nodes(fgrid) == num_elements_grid .+ (1,1)
    @test start(fgrid) == start_point
    @test finish(fgrid) == finish_point
    batch = rand(1:num_elements(fgrid)[2]); offset = rand(1:num_elements(fgrid)[1])
    cartesian_index_bench = batch * num_nodes(fgrid)[1] + offset
    cartesian_index_bench_to_test = CartesianIndex(offset , batch + 1)
    cartesian_to_test = cartesian_index(cartesian_index_bench , fgrid)
    @test cartesian_to_test == cartesian_index_bench_to_test

    # iterate over all fields of the object and check if are equal
    grid_to_test = grid(fgrid)
    fgrid_fields = fieldnames(typeof(grid_to_test))
    grid_test_bench, _ = _create_ferrite_rectangular_grid(
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
          (finish_point .+ start_point) ./ 2 ⊂ fgrid


    # test more specific ferrite functions
    Ωsets = _boundary_indexes(fgrid)

    # all sets must be unique
    @test length(unique(Ωsets)) == length(Ωsets)

    # test start point is in the first border cell
    start_cell = Set{FaceIndex}()
    push!(start_cell, FaceIndex((1, 1)))
    # TODO: Fix coordinates overlead
    # start_cell_coords = coordinates(start_cell, fgrid)
    # @test Vec{dimension(fgrid),node_type(fgrid)}(start_point) ∈ start_cell_coords

    # test corners are in border points
    points_per_segment = rand(2:5)
    Ωfgrid_points = border_points(fgrid, points_per_segment)
    bool_corners_in_borders = true
    for point in corners(fgrid)
        bool_corners_in_borders && point ∈ Ωfgrid_points
    end
    @test bool_corners_in_borders

    # test magnitude interpolation via ferrite grid
    magnitude = reshape(collect(1.0:prod(num_nodes(fgrid))), num_nodes(fgrid))
    fmagnitude_to_test = _convert_to_ferrite_nomenclature(magnitude, fgrid)
    fmagnitude_hand = [1, 2, 6, 5, 3, 7, 4, 8, 10, 9, 11, 12]
    @test fmagnitude_to_test == fmagnitude_hand

    # points that matches the grid nodes: start and final point value
    points_to_test = [start_point, finish_point]
    @test _interpolate(points_to_test, magnitude, fgrid) == [magnitude[1], magnitude[end]]
    # points that not matches the grid nodes:
    # consecutive in x + Δₓ / 2
    consecutive_x = start_point .+ (element_size(fgrid)[1], 0) ./ 2
    consecutive_x_hand = sum(magnitude[1,1] + magnitude[2,1]) / 2
    # consecutive in y + Δⱼ / 2
    consecutive_y = start_point .+ (0, element_size(fgrid)[2]) ./ 2
    consecutive_y_hand = sum(magnitude[1] + magnitude[1,2]) / 2
    # final point - (Δₓ, Δⱼ)/2
    final_cell_middle = finish_point .- element_size(fgrid)./2
    final_cell_middle_hand = sum(magnitude[end - 1:end,end - 1:end]) / 4

    points_to_test = [consecutive_x, consecutive_y, final_cell_middle]
    inter_to_test = _interpolate(points_to_test, magnitude, fgrid)
    inter_bench = [consecutive_x_hand, consecutive_y_hand, final_cell_middle_hand ]
    @test  inter_to_test == inter_bench




end
