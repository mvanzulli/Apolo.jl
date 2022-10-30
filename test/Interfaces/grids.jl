##############
# Grid tests #
##############

using Apolo.Geometry
using Apolo.Geometry: _convert_to_ferrite_nomenclature, _create_ferrite_rectangular_grid,
    _which_border, _corners, _borders, _boundary_indexes, _interpolate, _closest_point, _extrapolate

using Ferrite: Quadrilateral, Hexahedron, Grid, Set, FaceIndex, Vec
using Test: @test, @testset

const INTERVAL_START = LinRange(-10.0, 10.0, 20)
const INTERVAL_LENGTH = LinRange(1.0, 10.0, 20)
const TOLERANCE = 1e-3

@testset "Quadrilateral 2D Ferrite grids unitary tests" begin

    # Define a random grid
    start_point = (rand(INTERVAL_START), rand(INTERVAL_START))
    finish_point = start_point .+ (rand(INTERVAL_LENGTH), rand(INTERVAL_LENGTH))
    num_elements_grid = (3, 2)
    elemtype = Quadrilateral
    fgrid = FerriteStructuredGrid(start_point, finish_point, num_elements_grid, elemtype)


    # Check grid attributes with a handmade contrsuction
    grid_to_test = grid(fgrid)
    fgrid_fields = fieldnames(typeof(grid_to_test))
    grid_test_bench, _ = _create_ferrite_rectangular_grid(
        start_point,
        finish_point,
        num_elements_grid,
        elemtype
    )

    grids_are_the_same = true
    for field in fgrid_fields
        grids_are_the_same &&
            getproperty(grid_test_bench, field) == getproperty(grid_to_test, field)
    end
    @test grids_are_the_same


    # test getter functions
    @test corners(fgrid) == _corners(start_point, finish_point)
    @test extrema(fgrid) == ((start_point[1], finish_point[1]), (start_point[end], finish_point[end]))
    @test dimension(fgrid) == length(start_point)
    @test element_size(fgrid) == (finish_point .- start_point) ./ num_elements_grid
    @test element_type(fgrid) == elemtype
    @test node_type(fgrid) == eltype(finish_point)
    @test num_elements(fgrid) == num_elements_grid
    @test num_nodes(fgrid) == num_elements_grid .+ (1, 1)
    @test start(fgrid) == start_point
    @test finish(fgrid) == finish_point

    # test CartesianIndex computation
    cartindex_y = rand(1:num_nodes(fgrid)[2])
    cartindex_x = rand(1:num_nodes(fgrid)[1])
    cartesian_index_bench = (cartindex_y - 1) * num_nodes(fgrid)[1] +
                            cartindex_x
    cartesian_index_bench_to_test = CartesianIndex(cartindex_x, cartindex_y)
    cartesian_to_test = cartesian_index(cartesian_index_bench, fgrid)
    @test cartesian_to_test == cartesian_index_bench_to_test

    # test ∈
    @test start_point ⊂ fgrid && finish_point ⊂ fgrid &&
          (finish_point .+ start_point) ./ 2 ⊂ fgrid
    @test finish_point .+ (1.0, 1.0) ⊄ fgrid


    # test border functions with ferrite
    borders_bench = ["top", "bottom", "left", "right"]
    borders_to_test = _borders(fgrid)
    borders_are_the_same = true
    [borders_are_the_same && b in borders_bench for b in borders_to_test]
    @test borders_are_the_same

    # all sets must be unique
    Ωsets = _boundary_indexes(fgrid)
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
    fmagnitude_to_test, cartesian_to_test = _convert_to_ferrite_nomenclature(magnitude, fgrid)
    fmagnitude_hand = [1, 2, 6, 5, 3, 7, 4, 8, 10, 9, 11, 12]
    cartesian_indexes_hand = cartesian_index.(fmagnitude_hand, Ref(fgrid))
    @test fmagnitude_to_test == fmagnitude_hand
    @test cartesian_to_test == cartesian_indexes_hand

    # points that matches the grid nodes: start and final point value
    points_to_test = [start_point, finish_point]
    interpol_to_test = _interpolate(points_to_test, magnitude, :magnitude, fgrid)
    interpol_bench = [magnitude[1], magnitude[end]]
    @test interpol_to_test ≈ interpol_bench atol = TOLERANCE

    # points that not matches the grid nodes:
    # consecutive in x + Δₓ / 2
    consecutive_x = start_point .+ (element_size(fgrid)[1], 0) ./ 2
    consecutive_x_hand = sum(magnitude[1, 1] + magnitude[2, 1]) / 2
    # consecutive in y + Δⱼ / 2
    consecutive_y = start_point .+ (0, element_size(fgrid)[2]) ./ 2
    consecutive_y_hand = sum(magnitude[1] + magnitude[1, 2]) / 2
    # final point - (Δₓ, Δⱼ)/2
    final_cell_middle = finish_point .- element_size(fgrid) ./ 2
    final_cell_middle_hand = sum(magnitude[end-1:end, end-1:end]) / 4
    # test magnitude values
    points_to_test = [consecutive_x, consecutive_y, final_cell_middle]
    inter_to_test = _interpolate(points_to_test, magnitude, :magnitude, fgrid)
    inter_bench = [consecutive_x_hand, consecutive_y_hand, final_cell_middle_hand]
    @test inter_to_test ≈ inter_bench atol = TOLERANCE

    # test magnitude extrapolation via ferrite grid
    out_start = start_point .- (1.0, 1.0)
    out_finish = finish_point .+ (1.0, 1.0)
    out_consecutive_x = consecutive_x .- (0.0, 1.0)
    out_consecutive_y = consecutive_y .- (1.0, 0.0)
    @test out_start ⊄ fgrid
    @test out_finish ⊄ fgrid
    @test out_consecutive_x ⊄ fgrid
    @test out_consecutive_y ⊄ fgrid

    points_to_test = [out_start, out_finish, out_consecutive_x, out_consecutive_y]

    # test borders labels
    @test _which_border(out_consecutive_y, start_point, finish_point) == :left
    @test _which_border(out_consecutive_x, start_point, finish_point) == :bottom
    # test border closest points
    @test _closest_point(out_start, fgrid) == start_point
    @test _closest_point(out_finish, fgrid) == finish_point
    @test _closest_point(out_consecutive_x, fgrid) == (out_consecutive_x[1], start_point[2])
    @test _closest_point(out_consecutive_y, fgrid) == (start_point[1], out_consecutive_y[2])

    # test extrapolate
    extrapolation_to_test = _extrapolate(points_to_test, magnitude, fgrid)
    extrapolation_hand = [magnitude[1], magnitude[end], consecutive_x_hand, consecutive_y_hand]
    @test extrapolation_to_test ≈ extrapolation_hand atol = TOLERANCE
end

@testset "Hexahedron 3D Ferrite grids unitary tests" begin

    start_point = (rand(INTERVAL_START), rand(INTERVAL_START), rand(INTERVAL_START))
    start_point = (0.0, 0.0, 0.0)
    spacing_grid = (rand(INTERVAL_LENGTH), rand(INTERVAL_LENGTH), rand(INTERVAL_LENGTH))
    elements_size_grid = (1.0, 1.0, 1.0)
    spacing_grid = elements_size_grid
    num_elements_grid = (1, 3, 2)
    finish_point = start_point .+ spacing_grid .* elements_size_grid
    grid_dimension = length(spacing_grid)
    length_grid = finish_point .- start_point
    elemtype = Hexahedron


    fgrid = FerriteStructuredGrid(start_point, finish_point, num_elements_grid, elemtype)

    # Check grid attributes with a handmade contrsuction
    grid_to_test = grid(fgrid)
    fgrid_fields = fieldnames(typeof(grid_to_test))
    grid_test_bench, corners_bench = _create_ferrite_rectangular_grid(
        start_point,
        finish_point,
        num_elements_grid,
        elemtype
    )

    grids_are_the_same = true
    for field in fgrid_fields
        grids_are_the_same &&
            getproperty(grid_test_bench, field) == getproperty(grid_to_test, field)
    end
    @test grids_are_the_same

    # test getter functions
    @test corners(fgrid) == _corners(start_point, finish_point)
    @test extrema(fgrid) == (
        (start_point[1], finish_point[1]),
        (start_point[2], finish_point[2]),
        (start_point[3], finish_point[3])
    )
    @test dimension(fgrid) == length(start_point)
    @test num_elements(fgrid) == num_elements_grid
    @test num_nodes(fgrid) == num_elements_grid .+ (1, 1, 1)
    @test element_size(fgrid) == (finish_point .- start_point) ./ num_elements_grid
    @test element_type(fgrid) == elemtype
    @test node_type(fgrid) == eltype(finish_point)
    @test start(fgrid) == start_point
    @test finish(fgrid) == finish_point


    # test CartesianIndex computation
    cartindex_z = rand(1:num_nodes(fgrid)[3])
    cartindex_y = rand(1:num_nodes(fgrid)[2])
    cartindex_x = rand(1:num_nodes(fgrid)[1])
    cartesian_index_bench = (cartindex_z - 1) * prod(num_nodes(fgrid)[1:2]) +
                            (cartindex_y - 1) * num_nodes(fgrid)[1] +
                            cartindex_x
    cartesian_index_bench_to_test = CartesianIndex(cartindex_x, cartindex_y, cartindex_z)
    cartesian_to_test = cartesian_index(cartesian_index_bench, fgrid)
    @test cartesian_to_test == cartesian_index_bench_to_test

    # test ∈
    @test start_point ⊂ fgrid && finish_point ⊂ fgrid &&
          (finish_point .+ start_point) ./ 2 ⊂ fgrid
    @test (finish_point .+ (1.0, 1.0, 1.0)) ⊄ fgrid

    # test border functions with ferrite
    borders_bench = ["top", "bottom", "left", "right", "back", "front"]
    borders_to_test = _borders(fgrid)
    borders_are_the_same = true
    [borders_are_the_same && b in borders_bench for b in borders_to_test]
    @test borders_are_the_same

    # all sets must be unique
    Ωsets = _boundary_indexes(fgrid)
    @test length(unique(Ωsets)) == length(Ωsets)


    # test corners are in border points
    # TODO: EXTEND BORDER POINTS FOR 3D cases
    # points_per_segment = rand(2:5)
    # Ωfgrid_points = border_points(fgrid, points_per_segment)
    # bool_corners_in_borders = true
    # for point in corners(fgrid)
    # bool_corners_in_borders && point ∈ Ωfgrid_points
    # end
    # @test bool_corners_in_borders


    # test magnitude interpolation via ferrite grid
    magnitude = reshape(collect(1.0:prod(num_nodes(fgrid))), num_nodes(fgrid))
    fmagnitude_to_test, cartesian_to_test = _convert_to_ferrite_nomenclature(magnitude, fgrid)
    fmagnitude_hand = [1, 2, 4, 3, 9, 10, 12, 11, 6, 5, 14, 13, 8, 7, 16, 15, 17, 18, 20, 19, 22, 21, 24, 23]
    @test fmagnitude_to_test == fmagnitude_hand
    cartesian_indexes_hand = cartesian_index.(fmagnitude_hand, Ref(fgrid))
    @test fmagnitude_to_test == fmagnitude_hand
    @test cartesian_to_test == cartesian_indexes_hand


    # points that matches the grid nodes: start and final point value
    points_to_test = [start_point, finish_point]
    interpol_to_test = _interpolate(points_to_test, magnitude, :magnitude, fgrid)
    interpol_bench = [magnitude[1], magnitude[end]]
    @test interpol_to_test ≈ interpol_bench atol = TOLERANCE


    # points that not matches the grid nodes:
    # consecutive in x + Δₓ / 2
    consecutive_x = start_point .+ (element_size(fgrid)[1], 0, 0) ./ 2
    consecutive_x_hand = sum(magnitude[1, 1, 1] + magnitude[2, 1, 1]) / 2
    # consecutive in y + Δⱼ / 2
    consecutive_y = start_point .+ (0, element_size(fgrid)[2], 0) ./ 2
    consecutive_y_hand = sum(magnitude[1, 1, 1] + magnitude[1, 2, 1]) / 2
    # consecutive in y + Δₖ / 2
    consecutive_z = start_point .+ (0, 0, element_size(fgrid)[3]) ./ 2
    consecutive_z_hand = sum(magnitude[1, 1, 1] + magnitude[1, 1, 2]) / 2
    # final point - (Δₓ, Δⱼ)/2
    final_cell_middle = finish_point .- element_size(fgrid) ./ 2
    final_cell_middle_hand = sum(magnitude[end-1:end, end-1:end, end-1:end]) / 8

    # test magnitude values
    points_to_test = [consecutive_x, consecutive_y, consecutive_z, final_cell_middle]
    inter_to_test = _interpolate(points_to_test, magnitude, :magnitude, fgrid)
    inter_bench = [consecutive_x_hand, consecutive_y_hand, consecutive_z_hand, final_cell_middle_hand]
    @test inter_to_test ≈ inter_bench atol = TOLERANCE


    # test magnitude extrapolation via ferrite grid
    out_start = start_point .- (1.0, 1.0, 1.0)
    out_finish = finish_point .+ (1.0, 1.0, 1.0)
    out_consecutive_x = consecutive_x .- (0.0, 1.0, 0.0)
    out_consecutive_y = consecutive_y .- (1.0, 0.0, 0.0)
    out_consecutive_z = consecutive_z .- (1.0, 0.0, 0.0)
    out_top = (start_point[1] + (start_point[1] + finish_point[1]) / 2,
        finish_point[2] + 1.0,
        start_point[3] + (start_point[3] + finish_point[3]) / 2)
    out_bottom = (start_point[1] + (start_point[1] + finish_point[1]) / 2,
        start_point[2] - 1.0,
        start_point[3] + (start_point[3] + finish_point[3]) / 2)

    @test out_start ⊄ fgrid
    @test out_finish ⊄ fgrid
    @test out_top ⊄ fgrid
    @test out_bottom ⊄ fgrid
    @test out_consecutive_x ⊄ fgrid
    @test out_consecutive_y ⊄ fgrid
    @test out_consecutive_z ⊄ fgrid

    # test borders labels
    @test _which_border(out_consecutive_x, start_point, finish_point) == :bottom_front
    @test _which_border(out_consecutive_y, start_point, finish_point) == :left_front
    @test _which_border(out_consecutive_z, start_point, finish_point) == :left_bottom
    @test _which_border(out_top, start_point, finish_point) == :top
    @test _which_border(out_bottom, start_point, finish_point) == :bottom

    # test border closest points
    @test _closest_point(out_finish, fgrid) == finish_point
    @test _closest_point(out_start, fgrid) == start_point
    @test _closest_point(out_consecutive_x, fgrid) == (out_consecutive_x[1], start_point[2], start_point[3])
    @test _closest_point(out_consecutive_y, fgrid) == (start_point[1], out_consecutive_y[2], start_point[3])
    @test _closest_point(out_consecutive_z, fgrid) == (start_point[1], start_point[2], out_consecutive_z[3])
    @test _closest_point(out_top, fgrid) == (out_top[1], finish_point[2], out_top[3])
    @test _closest_point(out_bottom, fgrid) == (out_bottom[1], start_point[2], out_bottom[3])

    # test extrapolate
    points_to_test = [
        out_start, out_finish,
        out_consecutive_x, out_consecutive_y, out_consecutive_z
    ]
    extrapolation_to_test = _extrapolate(points_to_test, magnitude, fgrid)
    extrapolation_hand = [
        magnitude[1, 1, 1], magnitude[end, end, end],
        consecutive_x_hand, consecutive_y_hand, consecutive_z_hand
    ]
    @test extrapolation_to_test ≈ extrapolation_hand atol = TOLERANCE





end
