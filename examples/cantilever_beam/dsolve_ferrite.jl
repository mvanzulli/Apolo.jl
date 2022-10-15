# -----------------------------------
# Ferrite direct solver for example 2
# -----------------------------------
using BlockArrays, SparseArrays

struct LinearElasticityFerrite{T}
    G::T
    K::T
end

function dsolve_ferrite(linear_elasticity_prob, Ecand)

    # Extract the grid of fem data
    grid = linear_elasticity_prob.data_fem.grid
    # Extract the dof handler
    dh = linear_elasticity_prob.data_fem.dofs
    # Extract Dirichlet bcs
    dir = linear_elasticity_prob.data_fem.bcs.dirchlet

    # Create dirchlet boundary conditions
    function create_dbc!(dh, dirichlet_bc)
        dbc = Ferrite.ConstraintHandler(dh)
        Ferrite.add!(
            dbc,
            Ferrite.Dirichlet(
                dirichlet_bc.dof,
                dirichlet_bc.region,
                dirichlet_bc.imp_func,
                dirichlet_bc.imp_dofs,
            ),
        )
        Ferrite.close!(dbc)
        t = 0.0
        Ferrite.update!(dbc, t)
        return dbc
    end
    dbc = create_dbc!(
        linear_elasticity_prob.data_fem.dofs,
        linear_elasticity_prob.data_fem.bcs.dirchlet,
    )


    # Create cell values of dofs
    function create_values(interpolation_u, interpolation_p)
        # quadrature rules
        qr = Ferrite.QuadratureRule{2,Ferrite.RefTetrahedron}(3)
        face_qr = Ferrite.QuadratureRule{1,Ferrite.RefTetrahedron}(3)

        # geometric interpolation
        interpolation_geom = Ferrite.Lagrange{2,Ferrite.RefTetrahedron,1}()

        # cell and facevalues for u
        cellvalues_u = Ferrite.CellVectorValues(qr, interpolation_u, interpolation_geom)
        facevalues_u =
            Ferrite.FaceVectorValues(face_qr, interpolation_u, interpolation_geom)

        # cellvalues for p
        cellvalues_p = Ferrite.CellScalarValues(qr, interpolation_p, interpolation_geom)

        return cellvalues_u, cellvalues_p, facevalues_u
    end
    # define interpolations
    interpolation_u = Ferrite.Lagrange{2,Ferrite.RefTetrahedron,2}()
    interpolation_p = Ferrite.Lagrange{2,Ferrite.RefTetrahedron,1}()

    # crete values
    cellvalues_u, cellvalues_p, facevalues_u =
        create_values(interpolation_u, interpolation_p)


    # Assable local element and forces and add tensions
    function doassemble(
        cellvalues_u::Ferrite.CellVectorValues{dim},
        cellvalues_p::Ferrite.CellScalarValues{dim},
        facevalues_u::Ferrite.FaceVectorValues{dim},
        K::Ferrite.SparseMatrixCSC,
        grid::Ferrite.Grid,
        dh::Ferrite.DofHandler,
        mat,
    ) where {dim}

        f = zeros(Ferrite.ndofs(dh))
        assembler = Ferrite.start_assemble(K, f)
        nu = Ferrite.getnbasefunctions(cellvalues_u)
        np = Ferrite.getnbasefunctions(cellvalues_p)

        fe = BlockArrays.PseudoBlockArray(zeros(nu + np), [nu, np]) # local force vector
        ke = BlockArrays.PseudoBlockArray(zeros(nu + np, nu + np), [nu, np], [nu, np]) # local stiffness matrix

        # traction vector
        t = Ferrite.Vec{2}((0.0, -p))
        # cache ɛdev outside the element routine to avoid some unnecessary allocations
        ɛdev = [
            zero(Ferrite.SymmetricTensor{2,dim}) for
            i = 1:Ferrite.getnbasefunctions(cellvalues_u)
        ]

        for cell in Ferrite.CellIterator(dh)
            fill!(ke, 0)
            fill!(fe, 0)
            assemble_up!(
                ke,
                fe,
                cell,
                cellvalues_u,
                cellvalues_p,
                facevalues_u,
                grid,
                mp,
                ɛdev,
                t,
            )
            assemble!(assembler, celldofs(cell), fe, ke)
        end

        return K, f
    end


    function assemble_up!(
        Ke,
        fe,
        cell,
        cellvalues_u,
        cellvalues_p,
        facevalues_u,
        grid,
        mp,
        ɛdev,
        t,
    )

        n_basefuncs_u = Ferrite.getnbasefunctions(cellvalues_u)
        n_basefuncs_p = Ferrite.getnbasefunctions(cellvalues_p)
        u▄, p▄ = 1, 2
        Ferrite.reinit!(cellvalues_u, cell)
        Ferrite.reinit!(cellvalues_p, cell)

        # We only assemble lower half triangle of the stiffness matrix and then symmetrize it.
        @inbounds for q_point = 1:Ferrite.getnquadpoints(cellvalues_u)
            for i = 1:n_basefuncs_u
                ɛdev[i] = Ferrite.dev(
                    Ferrite.symmetric(Ferrite.shape_gradient(cellvalues_u, q_point, i)),
                )
            end
            dΩ = Ferrite.getdetJdV(cellvalues_u, q_point)
            for i = 1:n_basefuncs_u
                divδu = Ferrite.shape_divergence(cellvalues_u, q_point, i)
                δu = Ferrite.shape_value(cellvalues_u, q_point, i)
                for j = 1:i
                    Ke[BlockArrays.BlockIndex((u▄, u▄), (i, j))] +=
                        2 * mp.G * ɛdev[i] ⊡ ɛdev[j] * dΩ
                end
            end

            for i = 1:n_basefuncs_p
                δp = shape_value(cellvalues_p, q_point, i)
                for j = 1:n_basefuncs_u
                    divδu = shape_divergence(cellvalues_u, q_point, j)
                    Ke[BlockArrays.BlockIndex((p▄, u▄), (i, j))] += -δp * divδu * dΩ
                end
                for j = 1:i
                    p = shape_value(cellvalues_p, q_point, j)
                    Ke[BlockArrays.BlockIndex((p▄, p▄), (i, j))] += -1 / mp.K * δp * p * dΩ
                end

            end
        end

        symmetrize_lower!(Ke)

        # We integrate the Neumann boundary using the facevalues.
        # We loop over all the faces in the cell, then check if the face
        # is in our `"traction"` faceset.
        @inbounds for face = 1:nfaces(cell)
            if onboundary(cell, face) && (cellid(cell), face) ∈ getfaceset(grid, "traction")
                reinit!(facevalues_u, cell, face)
                for q_point = 1:getnquadpoints(facevalues_u)
                    dΓ = getdetJdV(facevalues_u, q_point)
                    for i = 1:n_basefuncs_u
                        δu = shape_value(facevalues_u, q_point, i)
                        fe[i] += (δu ⋅ t) * dΓ
                    end
                end
            end
        end
    end

    function symmetrize_lower!(K)
        for i = 1:size(K, 1)
            for j = i+1:size(K, 1)
                K[i, j] = K[j, i]
            end
        end
    end

    Emod = Ecand
    ν = 0.4
    Gmod = Emod / 2(1 + ν)
    Kmod = Emod * ν / ((1 + ν) * (1 - 2ν))
    mp = LinearElasticityFerrite(Gmod, Kmod)

    # assembly and solve
    K = Ferrite.create_sparsity_pattern(dh)
    K, f = doassemble(cellvalues_u, cellvalues_p, facevalues_u, K, grid, dh, mp)
    apply!(K, f, dbc)
    u = Symmetric(K) \ f
end
