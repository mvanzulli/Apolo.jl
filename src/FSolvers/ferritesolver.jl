
#####################################
# Main types for Ferrite.jl solver #
####################################

import Ferrite.DofHandler

using Apolo.Materials: SVK, value, lamé_params
using Apolo.ForwardProblem: ForwardProblemSolution, femdata, materials, dofsvals, label
import Apolo.ForwardProblem: _solve, _initialize!

using Reexport: @reexport
@reexport using Ferrite
using BlockArrays: BlockIndex, PseudoBlockArray
using SparseArrays: SparseMatrixCSC
using LinearAlgebra: Symmetric, norm

export FerriteForwardSolver


""" Struct with the interpolations for stress and displacements.

## Fields

- `itp_u` -- interpolation for displacement field u.
- `itp_σ` -- interpolation for stress σ field σ.
"""
struct InterStressDisp
    itp_u::Interpolation
    itp_σ::Interpolation
end

InterStressDisp(; itp_u=u_itp, itp_σ=σ_itp) = InterStressDisp(itp_u, itp_σ)

" Returns u interpolation. "
_u_itp(int::InterStressDisp) = int.itp_u

" Returns σ interpolation. "
_σ_itp(int::InterStressDisp) = int.itp_σ

" Struct with the quadrature rules for stress and displacements.

## Fields

- `face_u` -- quadrature rule for displacement field u across each face.
- `cell_u` -- quadrature rule for displacement field u inside each cell.
- `cell_σ` -- quadrature rule for stress field σ inside each cell.

"
struct QuadratureRulesStressDisp
    face_u::QuadratureRule
    cell_u::QuadratureRule
    cell_σ::QuadratureRule
end

"Constructor of QuadratureRulesStressDisp with heywords arguments"
function QuadratureRulesStressDisp(; face_u=face_u, cell_u=cell_u, cell_σ=cell_σ)
    return QuadratureRulesStressDisp(face_u, cell_u, cell_σ)
end

" Returns the quadrature rule for integration of the displacements field u inside a cell. "
_cell_u_qr(qrs::QuadratureRulesStressDisp) = qrs.cell_u

" Returns the quadrature rule for integration of the stress field σ inside a cell. "
_cell_σ_qr(qrs::QuadratureRulesStressDisp) = qrs.cell_σ

" Returns quadrature rules for integration of the displacements field u. "
_face_u_qr(qrs::QuadratureRulesStressDisp) = qrs.face_u

" Struct with the Face and Cell vector values for stress and displacements.

- `cell_val_u` -- cell vector value for u.
- `cell_val_σ` -- cell vector value for σ.

"
struct CellFaceValsStressDisp
    cell_val_u::Tuple{CellVectorValues,FaceVectorValues}
    cell_val_σ::CellScalarValues
end

"Constructor of QuadratureRulesStressDisp with heywords arguments"
CellFaceValsStressDisp(; cell_val_u=cvalu, cell_val_σ=cvalσ) = CellFaceValsStressDisp(cvalu, cvalσ)

" Retruns displacements cell values. "
_cell_values_u(cfvals::CellFaceValsStressDisp) = cfvals.cell_val_u[1]

" Retrun displacement face values. "
_face_values_u(cfvals::CellFaceValsStressDisp) = cfvals.cell_val_u[2]

" Retrun pressure face values. "
_cell_values_σ(cfvals::CellFaceValsStressDisp) = cfvals.cell_val_σ


" Ferrite solver struct.

## Fields

- `inter_geo `   -- geometric interpolations.
- `qrs`          -- quadrature rules.
- `inter_dofs`   -- interpolations defined for each dof.
- `dh`           -- ferrite `DofHandler` struct.
- `cellfacevals` -- cells face values.
- `nbasefuncs`   -- number of base functions

"
struct FerriteForwardSolver{IGEO,QRS,IDOFS,DH,CFV,NB} <: AbstractForwardProblemSolver
    inter_geo::IGEO
    qrs::QRS
    inter_dofs::IDOFS
    dh::DH
    cellfacevals::CFV
    nbasefuncs::NB

    "Constructor with user defined interpolations. "
    function FerriteForwardSolver(
        fproblem::AbstractForwardProblem,
        inter_geo::IGEO,
        qrs::QRS,
        inter_dofs::IDOFS,
    ) where {IGEO,QRS,IDOFS}

        # extract grid, dofs  and create: dof handler, cellface vals, nbasefuncs
        dofs_fproblem = dofs(fproblem)
        dh = DofHandler(grid(fproblem), dofs_fproblem, inter_dofs)
        cfv = create_values(dofs_fproblem, inter_dofs, inter_geo, qrs)
        nbasef = NumBaseFuncStressDisp(cfv)

        new{IGEO,QRS,IDOFS,typeof(dh),typeof(cfv),typeof(nbasef)}(
            inter_geo,
            qrs,
            inter_dofs,
            dh,
            cfv,
            nbasef,
        )
    end

    "`FerriteForwardSolver` constructor with default defined interpolations. "
    function FerriteForwardSolver(fproblem::AbstractForwardProblem)
        # load default elements and interpolations
        inter_geo, qrs, inter_dofs = default_elements(fproblem)
        # use the full inputs constructor
        FerriteForwardSolver(fproblem, inter_geo, qrs, inter_dofs)
    end

end

"Define default ferrite elements and interpolations"
function default_elements(fproblem)

    dimu = dimension(dofs(femdata(fproblem)).u)
    dimσ = dimension(dofs(femdata(fproblem)).σ)

    # elements to use
    tetra = RefTetrahedron

    # define interpolations
    inter_u = Lagrange{dimu,tetra,dimu}() # dim grid, element, order
    inter_σₓ = Lagrange{dimu,tetra,dimσ}()
    inter_dofs = InterStressDisp(itp_u=inter_u, itp_σ=inter_σₓ)

    # quadrature rules
    cell_qr = QuadratureRule{dimu,tetra}(3)
    face_qr = QuadratureRule{dimu - 1,tetra}(3)
    qrs = QuadratureRulesStressDisp(face_u=face_qr, cell_u=cell_qr, cell_σ=cell_qr)

    # geometric interpolation
    inter_geo = Lagrange{dimu,tetra,1}()

    return inter_geo, qrs, inter_dofs

end

# " Extract geometric interpolation of a ferrite forward problem solver. "
# getintgeo(ffs::FerriteForwardSolver) = ffs.inter_geo

# " Extract quadrature rules of a ferrite forward problem solver. "
# getqr(ffs::FerriteForwardSolver) = ffs.qrs

# " Extract quadrature rules of a ferrite forward problem solver. "
# getintdofs(ffs::FerriteForwardSolver) = ffs.inter_dofs

" Extract dofhandler."
dofhandler(ffs::FerriteForwardSolver) = ffs.dh

function dofhandler(sol::ForwardProblemSolution)
    if (:dh ∈ keys(sol.extra))
        sol.extra[:dh]
    else
        throw(ArgumentError("This solution type has no dofhandler"))
    end
end

" Extract cell face values."
cell_face_values(ffs::FerriteForwardSolver) = ffs.cellfacevals

# " Extract number of basis functions."
# getnbasesfuncs(ffs::FerriteForwardSolver) = ffs.nbasefuncs

" Creates a Ferrite dof handler "
function DofHandler(
    fgrid::FerriteStructuredGrid,
    dofs::StressDispDofs,
    inter_dofs::InterStressDisp,
)
    # Extract ferrite grid type and create dof handler
    dh = DofHandler(grid(fgrid))

    push!(dh, symbol(dofs.u), dimension(dofs.u), _u_itp(inter_dofs))
    push!(dh, symbol(dofs.σ), dimension(dofs.σ), _σ_itp(inter_dofs)) # check why is 1
    close!(dh)

    return dh

end


" Create cells and faces interpolation information for each dof "
function create_values(
    dofs,
    inter_dofs,
    inter_geo,
    qrs)

    # Cell and face values for u
    cellvalues_u = CellVectorValues(_cell_u_qr(qrs), _u_itp(inter_dofs), inter_geo)
    facevalues_u = FaceVectorValues(_face_u_qr(qrs), _u_itp(inter_dofs), inter_geo)

    # Cell values for p
    @assert dimension(dofs.σ) == 1 "TODO: Gerealize for vectorial σ dof"
    cellvalues_p = CellScalarValues(_cell_σ_qr(qrs), _σ_itp(inter_dofs), inter_geo)

    # Create struct
    CellFaceValsStressDisp((cellvalues_u, facevalues_u), cellvalues_p)
end

"Number of basis functions for displacements and tension"
struct NumBaseFuncStressDisp
    u::Int
    σ::Int
    function NumBaseFuncStressDisp(cfvals::CellFaceValsStressDisp)
        ndof_u = getnbasefunctions(_cell_values_u(cfvals))
        ndof_σ = getnbasefunctions(_cell_values_σ(cfvals))
        new(ndof_u, ndof_σ)
    end
end

"Initialize ferrite.jl forward problem solver"
function _initialize!(
    fproblem::FP,
    solver::FerriteForwardSolver,
    args...;
    kwargs...
) where {FP<:AbstractForwardProblem}

    # extract bcs and dofh
    bcs = boundary_conditions(fproblem)
    dh = dofhandler(solver)

    # create dirichlet_bcs
    dbc = create_dirichlet_bc(dh, bcs)

    # push dbc
    push!(fproblem.aux, :dbc => dbc)

    return nothing
end

""" Struct that contains the cell and face vals of a dof
### Fields:
- `u` -- cell and face values for displacements u
- `p` -- cell and face values for displacements u
"""


" Create Ferrite Dirichlet boundary conditions. "
function create_dirichlet_bc(
    dh::DofHandler,
    bcs::Dict{AbstractBoundaryCondition,Function},
)
    # create a constrain dof handler
    dbc = ConstraintHandler(dh)

    # iterate over bcs keys and get Dirichlet bc
    for bcᵢ in keys(bcs)
        if typeof(bcᵢ) == DirichletBC # dispatch or method to extract dirichlet bcs
            add!(
                dbc,
                Ferrite.Dirichlet(symbol(dofs(bcᵢ)),
                    getfaceset(dh.grid, string(label(bcᵢ))),
                    values_function(bcᵢ),
                    convert.(Int64, component(bcᵢ))),
            )
        end
    end
    # close dirichlet boundary condition dof handler
    close!(dbc)
    # update time
    t = 0.0
    update!(dbc, t)
    return dbc
end

"Solves a Linear Elasticty  problems with ferrite solver."
function _solve(
    fproblem::LinearElasticityProblem,
    solver::FerriteForwardSolver,
)

    # extract fproblem data
    ferrite_grid_fp = grid(grid(fproblem))
    mats = materials(fproblem)
    bcs = boundary_conditions(fproblem)
    cfv = cell_face_values(solver)
    cellvalues_u = _cell_values_u(cfv)
    facevalues_u = _face_values_u(cfv)
    cellvalues_σ = _cell_values_σ(cfv)

    dbc = fproblem.aux[:dbc]

    # unwarp solver
    dh = dofhandler(solver)
    K = create_sparsity_pattern(dh)

    # assamble the linear system
    K, f = doassemble(
        cellvalues_u,
        facevalues_u,
        cellvalues_σ,
        K,
        ferrite_grid_fp,
        dh,
        mats,
        bcs)

    # apply boundary condtions
    apply!(K, f, dbc)

    # solve u
    u = Symmetric(K) \ f

    # build solution
    sol = ForwardProblemSolution(
        solver,
        femdata(fproblem),
        mats,
        dofs(fproblem),
        u,
        Dict(:dh => dh),
    )

    return sol
end


" Starts Ferrite assembler and create sparsity pattern for u and p dof"
function doassemble(
    cellvalues_u::CellVectorValues{dim},
    facevalues_u::FaceVectorValues{dim},
    cellvalues_p::CellScalarValues{dim},
    K::SparseMatrixCSC,
    grid::Grid,
    dh::DofHandler,
    mats::Dict{AbstractMaterial,Function},
    bcs::Dict{AbstractBoundaryCondition,Function},
) where {dim}

    # external force for all dofs
    f = zeros(ndofs(dh))

    # start assemble process
    assembler = start_assemble(K, f)

    # get number of base functions for each dof
    nu = getnbasefunctions(cellvalues_u)
    np = getnbasefunctions(cellvalues_p)

    fe = PseudoBlockArray(zeros(nu + np), [nu, np]) # local force vector
    ke = PseudoBlockArray(zeros(nu + np, nu + np), [nu, np], [nu, np]) # local stiffness matrix

    # initialize in cache deformation tensor
    ɛdev = [zero(SymmetricTensor{2,dim}) for i = 1:getnbasefunctions(cellvalues_u)]

    # Iterate over each cell and assemble
    for cell in CellIterator(dh)
        fill!(ke, 0)
        fill!(fe, 0)

        # assemble element ke and f
        assemble_up!(
            ke,
            fe,
            cell,
            cellvalues_u,
            cellvalues_p,
            facevalues_u,
            grid,
            mats,
            ɛdev,
            bcs,
        )
        assemble!(assembler, celldofs(cell), fe, ke)
    end
    return K, f
end

" Assembles the tangent matrix and the force vector of the element "
function assemble_up!(
    Ke,
    fe,
    cell,
    cellvalues_u,
    cellvalues_p,
    facevalues_u,
    grid,
    mats,
    ɛdev,
    bcs,
)

    reinit!(cellvalues_u, cell)
    reinit!(cellvalues_p, cell)

    # find material of the cell
    matcell = material_cell(cell, mats, grid)

    # fill ke
    fill_linear_elasticKe!(
        Ke,
        matcell,
        ɛdev,
        cellvalues_u,
        cellvalues_p,
    )

    # fill fe
    applyNeumannBC!(fe, grid, bcs, cell, cellvalues_u, facevalues_u)

end

"Add Neumann boundary condition to the external force vector."
function applyNeumannBC!(
    fe,
    grid,
    bcs,
    cell,
    cellvalues_u,
    facevalues_u,
)
    # update cell values for each dof
    n_basefuncs_u = getnbasefunctions(cellvalues_u)

    # We integrate the Neumann boundary using the facevalues.
    # We loop over all the faces in the cell, then check if the face
    # Add NeumannLoadBC to the problem
    for bc in keys(bcs)

        if bc isa NeumannLoadBC
            # get label
            label_bc = string(label(bc))
            # iterate over each face cell
            @inbounds for face = 1:nfaces(cell)
                if onboundary(cell, face) && (cellid(cell), face) ∈ getfaceset(grid, label_bc)
                    reinit!(facevalues_u, cell, face)
                    for q_point = 1:getnquadpoints(facevalues_u)
                        dΓ = getdetJdV(facevalues_u, q_point)
                        for i = 1:n_basefuncs_u
                            # compute δu for the virtual work
                            δu = shape_value(facevalues_u, q_point, i)
                            # build tension vector
                            time = 0
                            t = Vec{length(bc.dir)}(bc.vals_func(time) .* Tuple(bc.dir)) # eval at time t = 0                            # compute nodal virtual work
                            fe[i] += (δu ⋅ t) * dΓ
                        end
                    end
                end
            end
        end
    end
end

"Fill linear elasitc tangent matrix of the element Ke."
function fill_linear_elasticKe!(
    Ke,
    matcell::SVK,
    ɛdev,
    cellvalues_u,
    cellvalues_p)

    # compute G and K
    Gmod, Kmod = lamé_params(matcell) #TODO move to materials interface from Hook to Bulk

    # update cell values for each dof
    n_basefuncs_u = getnbasefunctions(cellvalues_u)
    n_basefuncs_p = getnbasefunctions(cellvalues_p)
    u▄, p▄ = 1, 2

    # We only assemble lower half triangle of the stiffness matrix and then symmetrize it.
    @inbounds for q_point = 1:getnquadpoints(cellvalues_u)
        for i = 1:n_basefuncs_u
            ɛdev[i] = dev(symmetric(shape_gradient(cellvalues_u, q_point, i)))
        end
        dΩ = getdetJdV(cellvalues_u, q_point)
        for i = 1:n_basefuncs_u
            divδu = shape_divergence(cellvalues_u, q_point, i)
            δu = shape_value(cellvalues_u, q_point, i)
            for j = 1:i
                Ke[BlockIndex((u▄, u▄), (i, j))] += 2 * Gmod * ɛdev[i] ⊡ ɛdev[j] * dΩ
            end
        end

        for i = 1:n_basefuncs_p
            δp = shape_value(cellvalues_p, q_point, i)
            for j = 1:n_basefuncs_u
                divδu = shape_divergence(cellvalues_u, q_point, j)
                Ke[BlockIndex((p▄, u▄), (i, j))] += -δp * divδu * dΩ
            end
            for j = 1:i
                p = shape_value(cellvalues_p, q_point, j)
                Ke[BlockIndex((p▄, p▄), (i, j))] += -1 / Kmod * δp * p * dΩ
            end

        end
    end

    symmetrize_lower!(Ke)

    return nothing
end

" Extract the material of a cell"
function material_cell(cell, mats, grid)
    # find material of the cell
    idx_cell = cellid(cell)
    # check gird cell sets and extract first material label and then mat
    for (matlabel, matcellset) in getcellsets(grid)
        if idx_cell ∈ matcellset
            for mat in keys(mats)
                matlabel == label(mat) && return mat
            end
        end
    end
end

" Simetrize the tangent matrix"
function symmetrize_lower!(K)
    for i = 1:size(K, 1)
        for j = i+1:size(K, 1)
            K[i, j] = K[j, i]
        end
    end
end

#TODO: Add dispatch with the solver type in sol
"Gets the displacements values for a given vector of tuples and offset."
function _eval_displacements(
    sol::ForwardProblemSolution{S},
    vec_points::Vector{NTuple{D,T}},
    offset::NTuple{D,T}
) where {S<:FerriteForwardSolver,D,T}

    # Create a vector of Vec adding an offset
    vec_points = [Vec(vec .+ offset) for vec in vec_points]
    # Extract the dofs values according to ferrites nomenclature
    dofvals_sol = dofsvals(sol)

    # Create and eval the dof handler
    dh = dofhandler(sol)
    ferrite_grid = grid(grid(sol))
    ph = PointEvalHandler(ferrite_grid, vec_points)

    return Ferrite.get_point_values(ph, dh, dofvals_sol, symbol(sol.dofs.u))

end
