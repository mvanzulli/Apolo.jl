##########################################################
# Main types and function to compute volumetric matrices #
##########################################################

# Import dependencies to overlead
# external
import LazySets: Hyperrectangle
import Base: extrema, ∈, ∩, ∪,abs
# Add libraries to use
# internal
using Apolo.ForwardProblem: ForwardProblemSolution, getdofsvals

# External
import LinearAlgebra: Adjoint

using LazySets: LineSegment,
    vertices_list,
    isdisjoint

using Ferrite: AbstractCell, Grid, Vec,
    getfaceset,
    generate_grid

using Images: feature_transform, distance_transform


# Export interface functions and types
export SolutionBoundary, create_rectangular_grid,
getsides, get_boundary_points_def, get_boundary_ref,
binary_mat_zone, get_elements, abs, jaccard,
mesh_hyperrectangle, mean_hausdorff, geometric_simillarty, ⊂



#==================================================
Geometric interface with Ferrite.jl
==================================================#

# Constant variables or types
const VecOrTuple = Union{Vector,Tuple}
const DEFAULT_INTEGRATION_POINTS = 100
const DEFAULT_BOUNDARY_POINTS = 30


"""

Creates a rectangular grid.

## Fields
- `dimgrid`  -- grid dimension (max 3)
- `nelem`    -- Tuple with number of elements on each direction
- `start`    -- Start point
- `finish`   -- End point
- `elemType` -- Ferrite element type

"""
function create_rectangular_grid(
    dimgrid::Int,
    nelem::D,
    start::T,
    finish::T,
    elemType::Type{ET}
) where
{D<:VecOrTuple,T<:VecOrTuple,ET<:AbstractCell}
    # check dimensions match
    dimgrid = check_grid_dims(dimgrid, nelem, start, finish)
    dimgrid !==2 && throw(DimensionMismatch("By the moment is only dim = 2"))
    #TODO: extend the  implementation to generic dim grids


    corners = [
        Vec{dimgrid}(Tuple(start)),
        Vec{dimgrid}((finish[1], start[2])),
        Vec{dimgrid}(Tuple(finish)),
        Vec{dimgrid}((start[1], finish[2])),
    ]
    # create grid using a generic ferrite element type
    return generate_grid(elemType, nelem, corners)
end

" Check grid dimensions"
function check_grid_dims(dimgrid::Int, nelem::D, start::T, finish::T) where
{D<:VecOrTuple,T<:VecOrTuple}
    !(length(nelem) == length(start) == length(finish) == dimgrid) &&
        throw(DimensionMismatch("Check grid dimensions"))
    return length(start)
end

" Get border facesets"
function get_boundary_indexes(grid::Grid)
    Ωgrid_indexes = Set{FaceIndex}()
    # get faces indexes
    for border in ["top", "bottom", "left", "right"]
        # get faces index
        faces_index = getfaceset(grid, border)
        [push!(Ωgrid_indexes, face) for face in faces_index]
    end
    return Ωgrid_indexes
end

" Get grid point coordinates"
function Ferrite.getcoordinates(grid::Grid, Ωgrid_indexes::Set{FaceIndex})
    # Ωgrid_coords = Vector{Vec{getdim(grid), Float64}}()
    Ωgrid_coords = []
    for face in Ωgrid_indexes
        coords_face = getcoordinates(grid, face)
        for point in coords_face
            if (point ∉ Ωgrid_coords)
                push!(Ωgrid_coords, point)
            end
        end
    end
    return Ωgrid_coords
end

"Get grid extrema"
function extrema(grid::Grid)
    # extract boundary indexes
    boundary_indexes = get_boundary_indexes(grid)
    # extract coordinates
    boundary_coords = getcoordinates(grid, boundary_indexes)
    # extract xmax xmin ymax ymin
    (xₘᵢₙ, xₘₐₓ) = extrema(getindex.(boundary_coords, 1))
    (yₘᵢₙ, yₘₐₓ) = extrema(getindex.(boundary_coords, 2))
    # compute extrema
    return (xₘᵢₙ, yₘᵢₙ, xₘₐₓ, yₘₐₓ)
end

"Computes the grid border points"
function get_boundary_points_ref(grid::Grid, num_points_border::Int)

    # get max min
    (xₘᵢₙ, yₘᵢₙ, xₘₐₓ, yₘₐₓ) = extrema(grid)

    Ωleft_points = [Vec((xₘᵢₙ, y)) for y in range(yₘᵢₙ, yₘₐₓ, length=num_points_border)]
    Ωtop_points = [Vec((x, yₘₐₓ)) for x in range(xₘᵢₙ, xₘₐₓ, length=num_points_border)]
    Ωright_points = reverse([Vec((xₘₐₓ, y)) for y in range(yₘᵢₙ, yₘₐₓ, length=num_points_border)])
    Ωbottom_points = [Vec((x, yₘᵢₙ)) for x in range(xₘᵢₙ, xₘₐₓ, length=num_points_border)]


    # merge them into one single set
    merge_points(points, new_points) = [(p ∉ points && push!(points, p)) for p in new_points]

    # initialize and merge
    Ωgrid_points = Vector{Vec{getdim(grid),Float64}}()
    merge_points(Ωgrid_points, Ωleft_points)
    merge_points(Ωgrid_points, Ωtop_points)
    merge_points(Ωgrid_points, Ωright_points)
    merge_points(Ωgrid_points, Ωbottom_points)

    return Ωgrid_points
end

# TODO: Dispatch with ferrite solution
"Get the boudndary points of a forward problem solution obtained using ferrite.jl"
function get_boundary_points_def(
    sol::ForwardProblemSolution,
    num_points_border::Int,
)
    # extract grid
    grid = getgrid(sol)
    # extract boundary points of the grid in the reference configuration
    Ω_grid_ref = get_boundary_points_ref(grid, num_points_border)
    # extract u at the grid boundary
    dofs = getdofs(sol)
    dofu = dofs.u
    Ω_grid_u = get_dof_point_values(sol, Ω_grid_ref, dofu)
    # return boundary displacements
    return Ω_grid_u + Ω_grid_ref
end

"Solution set struct"
struct SolutionBoundary{dim,T<:Real} # es un polygon de LazySets
    elements::Vector{Vec{dim,T}}
end
# TODO: Para laburar con lazysets
#convert(::Type{Polygon}, ::SolutionBoundary)
#TODO: Add dispatch with 2D and 3D methods with solution boundary types
#
"Get solution set elements "
get_elements(sol_boundary::SolutionBoundary) = sol_boundary.elements

"Returns the number of sides that a solution set has"
nsides(sol_boundary::SolutionBoundary) = length(get_elements(sol_boundary))

" SolutionBoundary constructor with a solution"
function SolutionBoundary(
    sol::ForwardProblemSolution{FerriteForwardSolv},
    num_points_border=DEFAULT_BOUNDARY_POINTS
    )
    Ω_grid_def = get_boundary_points_def(sol, num_points_border)
    SolutionBoundary(Ω_grid_def)
end

"Returns the set sides

NOTE: elements of the solset must be in order
"
function getsides(sol_boundary::SolutionBoundary)
    # get elements
    elems = get_elements(sol_boundary)
    # initialize vector of sides
    allsides = Vector{Vector{typeof(get_elements(sol_boundary)[1])}}()

    for s in 1:nsides(sol_boundary)
        # TODO: ir hasta nsize -1 y al final del loop return
        # TODO: construir LineSegment allside un vector de LineSegment
        if s == nsides(sol_boundary)
            push!(allsides, [elems[s], elems[1]])
            return allsides
        else
            push!(allsides, [elems[s], elems[s+1]])
        end
    end
end

#==================================================
Geometric interface with LazySets.jl
==================================================#
#TODO: Si genero Vector de LineSegment
#=
box_approximation
box_approximation(::VPolygon)
Hyperrectangle
=#

"Lazy Sets Hyperrectangle constructor with a solution boundary "
function Hyperrectangle(solution_boundary::SolutionBoundary)
    # get borders
    (xₘᵢₙ, xₘₐₓ, yₘᵢₙ, yₘₐₓ) = get_maxmin(solution_boundary)
    # build HyperRectangle
    Hyperrectangle(low=[xₘᵢₙ, yₘᵢₙ], high=[xₘₐₓ, yₘₐₓ])
end

" Builds a Lazy Set Hyperrectangle that envelopes both solutions regions "
function Hyperrectangle(solution_boundary₁::SolutionBoundary, solution_boundary₂::SolutionBoundary)
    # compute the box enveloping both solutions
    (x1ₘᵢₙ, x1ₘₐₓ, y1ₘᵢₙ, y1ₘₐₓ) = get_maxmin(solution_boundary₁)
    (x2ₘᵢₙ, x2ₘₐₓ, y2ₘᵢₙ, y2ₘₐₓ) = get_maxmin(solution_boundary₂)

    # get the max min of each box
    (xₘᵢₙ, xₘₐₓ) = extrema([x1ₘᵢₙ, x1ₘₐₓ, x2ₘᵢₙ, x2ₘₐₓ])
    (yₘᵢₙ, yₘₐₓ) = extrema([y1ₘᵢₙ, y1ₘₐₓ, y2ₘᵢₙ, y2ₘₐₓ])

    # build box enveloping both solutions
    return Hyperrectangle(low=[xₘᵢₙ, yₘᵢₙ], high=[xₘₐₓ, yₘₐₓ])

end

"Builds a Lazy Set Hyperrectangle that envelopes the solution region "
function Hyperrectangle(sol::ForwardProblemSolution)
    # builds the solution boundary
    solution_boundary = SolutionBoundary(sol)
    return Hyperrectangle(solution_boundary)
end

"Builds a Lazy Set Hyperrectangle that envelopes both solutions "
function Hyperrectangle(
    sol₁::ForwardProblemSolution,
    sol₂::ForwardProblemSolution,
    num_points_per_border::Int=DEFAULT_BOUNDARY_POINTS,
)
    # build solution boundaries
    solution_boundary₁ = SolutionBoundary(sol₁, num_points_per_border)
    solution_boundary₂ = SolutionBoundary(sol₂, num_points_per_border)

    # build and return hyperrectangle containing both
    return Hyperrectangle(solution_boundary₁, solution_boundary₂)
end


"Checks if it is inside the box of the solution boundary"
function isinbox(p, solution_boundary::SolutionBoundary)
    # Compute the check using my own implementation
    # get vertexes
    (xₘᵢₙ, xₘₐₓ, yₘᵢₙ, yₘₐₓ) = get_maxmin(solution_boundary)
    isin = !(p[1] < xₘᵢₙ || p[1] > xₘₐₓ || p[2] < yₘᵢₙ || p[2] > yₘₐₓ)

    # Compute the check using LazySets
    # build hyper_rectangle
    hyper = Hyperrectangle(solution_boundary)
    # check lazy sets computes the same boolean
    if ((collect(p) ∈ hyper) == isin)
        return isin
    else
        throw(ArgumentError("The LazySets is different form this pkg check"))
    end
end

"Get borders of a solution boundary set"
function get_maxmin(solution_boundary::SolutionBoundary)
    # get x, y points
    x_points = getindex.(get_elements(solution_boundary), 1)
    y_points = getindex.(get_elements(solution_boundary), 2)
    # get max min in x and y
    (xₘᵢₙ, xₘₐₓ) = extrema(x_points)
    (yₘᵢₙ, yₘₐₓ) = extrema(y_points)
    # return a named Tuple with vertex
    return (xₘᵢₙ, xₘₐₓ, yₘᵢₙ, yₘₐₓ)
end

"Checks if a point is inside a solution boundary set using ray casting "
function ray_casting(p, solution_boundary::SolutionBoundary)
    # Sketch:
    # --------
    # |  x(p)|
    # |      |       x(pout)
    # -------

    # select a point outside the solution set
    (xₘᵢₙ, xₘₐₓ, yₘᵢₙ, yₘₐₓ) = get_maxmin(solution_boundary)
    pout = (xₘₐₓ, yₘₐₓ) .+ (1.0, 1.0)
    # create a segment between p and pout
    segment_p_pout = LineSegment(collect(p), collect(pout))

    # compute sides of a solution set
    solset_sides = getsides(solution_boundary)

    intersections = 0
    # check for all sides
    for side in solset_sides
        # build segment
        side_segment = LineSegment(collect(side[1]), collect(side[2]))
        # intersections
        if (!isdisjoint(side_segment, segment_p_pout))
            intersections += 1
        end
    end

    return intersections & 1 == 1
end

"Checks if a point is inside a SolutionBoundary"
function ∈(p, solution_boundary::SolutionBoundary)
    # check if it inside the box and if not break
    !isinbox(p, solution_boundary) && return false
    # if is in the box then check using ray tracing
    ray_casting(p, solution_boundary::SolutionBoundary)
end

"Given two vectors xᵥ, yᵥ builds an image with 1 and 0 where a solution_boundary is contained"
function binary_mat_zone(xᵥ, yᵥ, solution_boundary::SolutionBoundary)
    # get xvec and yvec sizes
    (m, n) = (length(xᵥ), length(yᵥ))
    # initialize the matrix
    mat = Matrix{Bool}(undef, (m, n))
    for (index_y, y) in enumerate(yᵥ)
        for (index_x, x) in enumerate(xᵥ)
            mat[index_x, index_y] = ((x, y) ∈ solution_boundary)
        end
    end
    return mat
end

" Computes the binary mat zone given a specific hyper rectangle domain"
function binary_mat_zone(
    solution_boundary::SolutionBoundary,
    hyper_rectangle::Hyperrectangle,
    num_inter_per_axis::Int=DEFAULT_INTEGRATION_POINTS,
    )

    # mesh the hyper rectangle with `num_inter_per_axis` per axis
    (xᵥ, yᵥ, Δx, Δy) = mesh_hyperrectangle(hyper_rectangle, num_inter_per_axis)

    # build the binary matrix
    bin_mat = binary_mat_zone(xᵥ, yᵥ, solution_boundary)
    return bin_mat, Δx, Δy
end

"Create a structured two x and y vector to integrate inside a box "
function mesh_hyperrectangle(
    hyperrect::Hyperrectangle,
    num_inter_per_axis::Int,
)
    # computes the hyper extrema
    ((xₘᵢₙ, yₘᵢₙ), (xₘₐₓ, yₘₐₓ)) = extrema(hyperrect)

    # create the grid where both solutions are integrated
    # in order to use EDV function the envelope box must comply Δx/Δy = 1

    # Build x vector
    num_inter_per_axis_x = num_inter_per_axis
    Δx = (xₘₐₓ - xₘᵢₙ) / num_inter_per_axis
    xᵥ = range(xₘᵢₙ + Δx/2 , xₘₐₓ - Δx/2 , length = num_inter_per_axis_x)

    # Build y vector
    # TODO: mesh with different number of pixels in x and y
    # select Δy
    Δy = Δx
    # compute the number of point to cover yₘᵢₙ
    num_inter_per_axis_y = cld( (yₘₐₓ - Δy/2) - (yₘᵢₙ + Δy/2), Δy)
    # compute the number of point to cover yₘᵢₙ
    yₘₐₓ_mesh = (yₘᵢₙ + Δy/2) + num_inter_per_axis_y * Δy
    # build the
    yᵥ =  (yₘᵢₙ + Δy/2) : Δy : yₘₐₓ_mesh

    return (xᵥ, yᵥ, Δx, Δy)
end

# " Computes two solution boundary intersection matrix (1 ∈ both 0 if not) and area"
function ∩(
    solution_boundary₁::SolutionBoundary,
    solution_boundary₂::SolutionBoundary,
    num_inter_per_axis::Int=DEFAULT_INTEGRATION_POINTS,
)

    # build the hyperrectangle containing both solutions
    hyperect_solutions = Hyperrectangle(solution_boundary₁, solution_boundary₂)

    # mesh the hyper rectangle with `num_inter_per_axis` per axis
    (xᵥ, yᵥ, Δx, Δy) = mesh_hyperrectangle(hyperect_solutions, num_inter_per_axis)

    # compute the intersection matrix
    intersection_mat = binary_mat_zone(xᵥ, yᵥ, solution_boundary₁) .* binary_mat_zone(xᵥ, yᵥ, solution_boundary₂)

    # compute the ∩ area integrating the intersection_mat boolean matrix inside
    # the box
    intersection_area = ∫(intersection_mat, Δx, Δy)

    return intersection_mat, intersection_area
end

"Computes the integral of boolean matrix in a domain with a [Δx x Δy] voxel "
∫(mat, Δx, Δy) = sum(mat) * Δx * Δy

" Computes two solution boundary intersection matrix (1 ∈ both 0 if not)"
function ∩(
    sol₁::ForwardProblemSolution,
    sol₂::ForwardProblemSolution,
    num_points_per_border::Int=DEFAULT_BOUNDARY_POINTS,
    num_inter_per_axis::Int=DEFAULT_INTEGRATION_POINTS,
)

    # build solution boundaries
    solution_boundary₁ = SolutionBoundary(sol₁, num_points_per_border)
    solution_boundary₂ = SolutionBoundary(sol₂, num_points_per_border)

    # compute ∩ matrix
    return ∩(solution_boundary₁, solution_boundary₂, num_inter_per_axis)
end

"Computes the solution boundary area or |solution_boundary|"
function abs(
    solution_boundary::SolutionBoundary,
    num_inter_per_axis::Int=DEFAULT_INTEGRATION_POINTS
)

    # build solution_boundary hyperrectangle
    hyperrectangle_solution = Hyperrectangle(solution_boundary)

    # mesh hyperrectangle solution
    (xᵥ, yᵥ, Δx, Δy) = mesh_hyperrectangle(hyperrectangle_solution, num_inter_per_axis)

    # compute the intersection matrix
    intersection_mat = binary_mat_zone(xᵥ, yᵥ, solution_boundary)

    # compute the area
    return ∫(intersection_mat, Δx, Δy)

end

"Computes the solution boundary area or |solution_boundary|"
function abs(
    sol::ForwardProblemSolution,
    num_points_per_border::Int=DEFAULT_BOUNDARY_POINTS,
    num_inter_per_axis::Int=DEFAULT_INTEGRATION_POINTS
)

    # build solution boundary
    solution_boundary = SolutionBoundary(sol, num_points_per_border)

    # call abs for a solution boundary
    return abs(solution_boundary, num_inter_per_axis)

end

"Returns the union beteween two different solutions"
function ∪(
    sol₁::ForwardProblemSolution,
    sol₂::ForwardProblemSolution,
    num_points_per_border::Int=DEFAULT_BOUNDARY_POINTS,
)

    # abs of each solution
    abs_sol₁ = abs(sol₁, num_points_per_border, num_inter_per_axis)
    abs_sol₂ = abs(sol₂, num_points_per_border, num_inter_per_axis)

    # computes intersection area
    _,intersection_area = ∩(sol₁, sol₂, num_points_per_border, num_inter_per_axis)

    # return the union result
    return ((abs_sol₁ + abs_sol₂) - intersection_area )
end

"Returns Jaccard's similarity number between two different solutions"
function jaccard(
    sol₁::ForwardProblemSolution,
    sol₂::ForwardProblemSolution,
    num_points_per_border::Int=DEFAULT_BOUNDARY_POINTS,
    num_inter_per_axis::Int=DEFAULT_INTEGRATION_POINTS,
)

    # computes intersection area
    _,intersection_area = ∩(sol₁, sol₂, num_points_per_border, num_inter_per_axis)

    # abs of each solution
    abs_sol₁ = abs(sol₁, num_points_per_border, num_inter_per_axis)
    abs_sol₂ = abs(sol₂, num_points_per_border, num_inter_per_axis)

    # return Jaccard's number
    return intersection_area / ((abs_sol₁ + abs_sol₂) - intersection_area )
end

"Returns Dice similarity number between two different solutions"
function dice(
    sol₁::ForwardProblemSolution,
    sol₂::ForwardProblemSolution,
    num_points_per_border::Int=DEFAULT_BOUNDARY_POINTS,
    num_inter_per_axis::Int=DEFAULT_INTEGRATION_POINTS,
)
    # computes intersection area
    _,intersection_area = ∩(sol₁, sol₂, num_points_per_border, num_inter_per_axis)

    # abs of each solution
    abs_sol₁ = abs(sol₁, num_points_per_border, num_inter_per_axis)
    abs_sol₂ = abs(sol₂, num_points_per_border, num_inter_per_axis)

    # return Jaccard's number
    return 2 * intersection_area / (abs_sol₁ + abs_sol₂)
end



"Returns the Hausdorff_distance between two forward problem solutions"
function mean_hausdorff(
    sol₁::ForwardProblemSolution,
    sol₂::ForwardProblemSolution,
    num_points_per_border::Int=DEFAULT_BOUNDARY_POINTS,
    num_inter_per_axis::Int=DEFAULT_INTEGRATION_POINTS,
    )

    # build the hyperrectangle containing both solutions
    hyperect_solutions = Hyperrectangle(sol₁, sol₂)

    # build solutions boundary sets
    solution_boundary₁ = SolutionBoundary(sol₁, num_points_per_border)
    solution_boundary₂ = SolutionBoundary(sol₂, num_points_per_border)

    # intersection_matrices
    intersection_mat₁, Δx₁, Δy₁ = binary_mat_zone(
        solution_boundary₁,
        hyperect_solutions,
        num_inter_per_axis,
        )

    intersection_mat₂, Δx₂, Δy₂ = binary_mat_zone(
        solution_boundary₂,
        hyperect_solutions,
        num_inter_per_axis,
        )

    @assert Δx₁ == Δy₁ == Δx₂ == Δy₂

    # Build cartesian index matrix and EDV
    # solution 1:
    cart_index_matrix₁ = feature_transform(intersection_mat₁)
    # euclidean distance variation (EDV)
    EDV₁ = Δx₁ * distance_transform(cart_index_matrix₁) #TODO: Fix me if Δx / Δy is not 1
    # solution 2:
    cart_index_matrix₂ = feature_transform(intersection_mat₂)
    # euclidean distance variation (EDV)
    EDV₂ = Δx₂ * distance_transform(cart_index_matrix₂)


    # Find indexes of each solution in the image ([i,j] such that EDV[i,j] == 0)
    # solution 1:
    sol₁_cart_index = findall(EDV₁ .== 0)
    # number of pixels inside solution 1
    N₁ = length(sol₁_cart_index)
    # solution 2:
    sol₂_cart_index = findall(EDV₂ .== 0)
    # number of pixels inside solution 2
    N₂ = length(sol₂_cart_index)

    # hausdorff distances
    hsol₂ =  1/N₂ * sum(EDV₂[sol₁_cart_index])
    hsol₁ =  1/N₁ * sum(EDV₁[sol₂_cart_index])



    return 1/2 * (hsol₁ + hsol₂)

end
""

"Computes geometric similarity based on S.Martinez and M. J. Rpuerez 2018
between two forward problems"
function geometric_simillarty(
    sol₁::ForwardProblemSolution,
    sol₂::ForwardProblemSolution,
    num_points_per_border::Int=DEFAULT_BOUNDARY_POINTS,
    num_inter_per_axis::Int=DEFAULT_INTEGRATION_POINTS,
)

    jc = jaccard(sol₁, sol₂, num_points_per_border, num_inter_per_axis)
    mhd = mean_hausdorff(sol₁, sol₂, num_points_per_border, num_inter_per_axis)

    return log((1 - jc) *  mhd)
end

#TODO: add compute metrix with only one build of silution boundary struct
