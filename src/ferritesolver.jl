
#####################################
# Main types for Ferrite.jl solver #
####################################

# Import dependencies to overlead
import IdenGPU.ForwardProblem: _solve, initialize!

# Add libraries to use
# Internal
using IdenGPU.Materials: SVK
using IdenGPU.ForwardProblem: ForwardProblemSolution, getfemdata, getmats, getdofsvals
# External
using BlockArrays: BlockIndex, PseudoBlockArray
using SparseArrays: SparseMatrixCSC
using LinearAlgebra: Symmetric
using Ferrite

# Export interface functions
export InterStressDisp, getuinter, getσinter,
    QuadratureRulesStressDisp, getcellqr, getfaceqr,
    CellFaceValsStressDisp, getucellval, getufaceval, getσcellval,
    FerriteForwardSolv, getdh, getintgeo, getintdofs, create_dirichlet_bc,
    create_values, get_dof_point_values

" Struct with the interpolations for stress and displacements."
Base.@kwdef struct InterStressDisp
    u::Interpolation
    σ::Interpolation
end

" get u interpolation. "
getuinter(ints::InterStressDisp) = ints.u

" get σ interpolation. "
getσinter(ints::InterStressDisp) = ints.σ

" Struct with the quadrature rules for stress and displacements." # Is assumed that the same qrule is used for face cell
Base.@kwdef struct QuadratureRulesStressDisp
    faceu::QuadratureRule
    cellu::QuadratureRule
    cellσ::QuadratureRule
end

" Get cell u quadrature rule. "
getcelluqr(qrs::QuadratureRulesStressDisp) = qrs.cellu
" Get face u quadrature rule. "
getfaceuqr(qrs::QuadratureRulesStressDisp) = qrs.faceu
" Get cell σ quadrature rule. "
getcellσqr(qrs::QuadratureRulesStressDisp) = qrs.cellq

" Struct with the Face and Cell vector values for stress and displacements."
Base.@kwdef struct CellFaceValsStressDisp
    u::Tuple{CellVectorValues,FaceVectorValues}
    σ::CellScalarValues
end

" Retrun displacements cell values. "
getucellval(cfvals::CellFaceValsStressDisp) = cfvals.u[1]
" Retrun displacement face values. "
getufaceval(cfvals::CellFaceValsStressDisp) = cfvals.u[2]
" Retrun pressure face values. "
getσcellval(cfvals::CellFaceValsStressDisp) = cfvals.σ


" Ferrite solver struct. "
struct FerriteForwardSolv{IGEO,QRS,IDOFS,DH,CFV,NB} <: AbstractForwardProbSolver
    inter_geo::IGEO
    qrs::QRS
    inter_dofs::IDOFS
    dh::DH
    cellfacevals::CFV
    nbasefuncs::NB

    "`FerriteForwardSolv` constructor with user defined interpolations. "
    function FerriteForwardSolv(
        fproblem::AbstractForwardProblem,
        inter_geo::IGEO,
        qrs::QRS,
        inter_dofs::IDOFS,
    ) where {IGEO,QRS,IDOFS}
        # extract grid, dofs  and create: dof handler, cellface vals, nbasefuncs
        grid = getgrid(fproblem)
        dofs = getdofs(fproblem)
        dh = create_dofhandler(grid, dofs, inter_dofs)
        cfv = create_values(dofs, inter_dofs, inter_geo, qrs)
        nbasef = NumBaseFuncStressDisp(cfv)
        # construct
        new{IGEO,QRS,IDOFS,typeof(dh),typeof(cfv),typeof(nbasef)}(
            inter_geo,
            qrs,
            inter_dofs,
            dh,
            cfv,
            nbasef,
        )
    end
    "`FerriteForwardSolv` constructor with default defined interpolations. "
    function FerriteForwardSolv(fproblem::AbstractForwardProblem)
        # load default elements and interpolations
        inter_geo, qrs, inter_dofs = default_elements(fproblem)
        # use the full inputs constructor 
        FerriteForwardSolv(fproblem, inter_geo, qrs, inter_dofs)
    end

end

"Define default ferrite elements and interpolations"
function default_elements(fproblem)


    # Exract u, σ and grid dimensions
    dimu = getdim(fproblem.data.dofs.u)
    dimσ = getdim(fproblem.data.dofs.σ)
    dimgrid = getdim(fproblem.data.grid)

    # elements to use 
    tetra = RefTetrahedron
    # define interpolations
    inter_u = Lagrange{dimu,tetra,dimu}() # dim grid, element, order
    inter_σₓ = Lagrange{dimu,tetra,dimσ}()
    inter_dofs = InterStressDisp(u=inter_u, σ=inter_σₓ)

    # quadrature rules
    cell_qr = QuadratureRule{dimu,tetra}(3)
    face_qr = QuadratureRule{dimu - 1,tetra}(3)
    qrs = QuadratureRulesStressDisp(faceu=face_qr, cellu=cell_qr, cellσ=cell_qr)
    # geometric interpolation
    inter_geo = Lagrange{dimu,tetra,1}()

    return inter_geo, qrs, inter_dofs

end

" Extract geometric interpolation of a ferrite forward problem solver. "
getintgeo(ffs::FerriteForwardSolv) = ffs.inter_geo

" Extract quadrature rules of a ferrite forward problem solver. "
getqr(ffs::FerriteForwardSolv) = ffs.qrs

" Extract quadrature rules of a ferrite forward problem solver. "
getintdofs(ffs::FerriteForwardSolv) = ffs.inter_dofs

" Extract dofhandler."
getdh(ffs::FerriteForwardSolv) = ffs.dh

function getdh(sol::ForwardProblemSolution)
    if (:dh  ∈ keys(sol.extra) )
        sol.extra[:dh]
    else 
        throw(ArgumentError("This solution type has no dofhandler"))
    end
end
    
" Extract cell face values."
getcellfacevalues(ffs::FerriteForwardSolv) = ffs.cellfacevals

" Extract number of basis functions."
getnbasesfuncs(ffs::FerriteForwardSolv) = ffs.nbasefuncs

" Creates a Ferrite dof handler "
function create_dofhandler(
    grid::Grid,
    dofs::StressDispDofs,
    inter_dofs::InterStressDisp,
)
    # create a dof handle
    dh = DofHandler(grid)
    # Push displacement into dof handler
    push!(dh, getsym(dofs.u), getdim(dofs.u), inter_dofs.u)
    # Push pressure into dof handler
    push!(dh, getsym(dofs.σ), getdim(dofs.σ), inter_dofs.σ) # check why is 1 
    # Close dofhandler and return it
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
    cellvalues_u = CellVectorValues(qrs.cellu, inter_dofs.u, inter_geo)
    facevalues_u = FaceVectorValues(qrs.faceu, inter_dofs.u, inter_geo)

    # Cell values for p
    @assert getdim(dofs.σ) == 1 "TODO: Gerealize for vectorial σ dof"
    cellvalues_p = CellScalarValues(qrs.cellσ, inter_dofs.σ, inter_geo)

    # Create struct
    CellFaceValsStressDisp((cellvalues_u, facevalues_u), cellvalues_p)
end

"Number of basis functions for displacements and tension"
struct NumBaseFuncStressDisp
    u::Int
    σ::Int
    function NumBaseFuncStressDisp(cfvals::CellFaceValsStressDisp)
        ndof_u = getnbasefunctions(getucellval(cfvals))
        ndof_σ = getnbasefunctions(getσcellval(cfvals))
        new(ndof_u, ndof_σ)
    end
end

"Initialize ferrite.jl forward problem solver"
function initialize!(
    fproblem::FP,
    solver::FerriteForwardSolv,
    args...;
    kwargs...
) where {FP<:AbstractForwardProblem}

    # extract bcs and dofh
    bcs = getbcs(fproblem)
    dh = getdh(solver)

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
                Ferrite.Dirichlet(getsym(getdofs(bcᵢ)),
                    getfaceset(dh.grid, string(getlabel(bcᵢ))),
                    vals_func(bcᵢ),
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

"Solves linear elasiticy problems with ferrite solver."
function _solve(
    fproblem::LinearElasticityProblem,
    solver::FerriteForwardSolv,
)

    # unwarp fproblem data 
    grid = getgrid(fproblem)
    mats = getmats(fproblem)
    bcs = getbcs(fproblem)
    cfv = getcellfacevalues(solver)
    cellvalues_u = getucellval(cfv)
    facevalues_u = getufaceval(cfv)
    cellvalues_σ = getσcellval(cfv)

    dbc = fproblem.aux[:dbc]

    # unwarp solver
    dh = getdh(solver)
    K = create_sparsity_pattern(dh)

    # assamble the linear system
    K, f = doassemble(
        cellvalues_u,
        facevalues_u,
        cellvalues_σ,
        K,
        grid,
        dh,
        mats,
        bcs)

    # apply boundary condtions
    apply!(K, f, dbc)

    # solve u
    u = Symmetric(K) \ f

    # build solution
    sol = ForwardProblemSolution(solver, getfemdata(fproblem), mats, getdofs(fproblem), u, Dict(:dh => dh))

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
end;

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
    matcell = getcellmat(cell, mats, grid)
    (Gmod, Kmod) = getmatparams(matcell)

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

        if typeof(bc) == NeumannLoadBC
            # get label
            label = string(getlabel(bc))
            # iterate over each face cell
            @inbounds for face = 1:nfaces(cell)
                if onboundary(cell, face) && (cellid(cell), face) ∈ getfaceset(grid, label)
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
    (Gmod, Kmod) = getmatparams(matcell) #TODO move to materials interface from Hook to Bulk

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
function getcellmat(cell, mats, grid)
    # find material of the cell
    idx_cell = cellid(cell)
    # check gird cell sets and extract first material label and then mat 
    for (matlabel, matcellset) in getcellsets(grid)
        if idx_cell ∈ matcellset
            for mat in keys(mats)
                matlabel == getlabel(mat) && return mat
            end
        end
    end
end

" Extract svk material parameters to use with ferrite nomenclature"
function getmatparams(svk::SVK)
    # Extract E and ν parameters
    E = getval(svk[:E])
    ν = getval(svk[:ν])

    # Compute Lamé parameters 
    Gmod = E / 2(1 + ν)
    Kmod = E * ν / ((1 + ν) * (1 - 2ν))

    return Gmod, Kmod
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
" Get the values  "
function get_dof_point_values(sol::ForwardProblemSolution, vec_points::Vector{Vec{dim,T}}, dof) where {dim,T}
    # get all dofs u value
    dofvals = getdofsvals(sol)
    # get dof handler
    dh = getdh(sol)
    # create a point handler
    grid = getgrid(sol)
    ph = PointEvalHandler(grid, vec_points)
    # eval the points
   return Ferrite.get_point_values(ph, dh, dofvals, getsym(dof))

end




