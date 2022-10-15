# -------------------------------------------------------------------------------------
# This is the Ferrite script of a one dimensional identification problem of Zerpa 2019 
# -------------------------------------------------------------------------------------
using Ferrite, Infiltrator
using BlockArrays, SparseArrays, LinearAlgebra
# 
import Plots
using LaTeXStrings


Lᵢₛ = 2.0;
Lⱼₛ = 1.0;
Eᵣ = 14e6;
ν = 0.4;
p = 1e3;
"""
function to create the grid

"""
function create_CMAME1_grid(nx, ny, Lᵢₛ=Lᵢₛ, Lⱼₛ=Lⱼₛ)
    # set corners of the mesh
    corners =
        [Vec{2}((0.0, 0.0)), Vec{2}((Lᵢₛ, 0.0)), Vec{2}((Lᵢₛ, Lⱼₛ)), Vec{2}((0.0, Lⱼₛ))]
    # create grid using Quadrilateral  type of element
    grid = generate_grid(Triangle, (nx, ny), corners)
    # facesets for boundary conditions
    addfaceset!(grid, "clamped", x -> norm(x[1]) ≈ 0.0)
    addfaceset!(grid, "traction", x -> norm(x[1]) ≈ Lᵢₛ)
    return grid
end;

function create_values(interpolation_u, interpolation_p)
    # quadrature rules
    qr = QuadratureRule{2,RefTetrahedron}(3)
    face_qr = QuadratureRule{1,RefTetrahedron}(3)

    # geometric interpolation
    interpolation_geom = Lagrange{2,RefTetrahedron,1}()

    # cell and facevalues for u
    cellvalues_u = CellVectorValues(qr, interpolation_u, interpolation_geom)
    facevalues_u = FaceVectorValues(face_qr, interpolation_u, interpolation_geom)

    # cellvalues for p
    cellvalues_p = CellScalarValues(qr, interpolation_p, interpolation_geom)

    return cellvalues_u, cellvalues_p, facevalues_u
end;

function create_dofhandler(grid, ipu, ipp)
    dh = DofHandler(grid)
    push!(dh, :u, 2, ipu) # displacement
    push!(dh, :p, 1, ipp) # pressure
    close!(dh)
    return dh
end;

function create_bc(dh)
    dbc = ConstraintHandler(dh)
    Ferrite.add!(
        dbc,
        Dirichlet(:u, getfaceset(dh.grid, "clamped"), (x, t) -> zero(Vec{2}), [1, 2]),
    )
    close!(dbc)
    t = 0.0
    update!(dbc, t)
    return dbc
end;

struct LinearElasticity{T}
    G::T
    K::T
end

function doassemble(
    cellvalues_u::CellVectorValues{dim},
    cellvalues_p::CellScalarValues{dim},
    facevalues_u::FaceVectorValues{dim},
    K::SparseMatrixCSC,
    grid::Grid,
    dh::DofHandler,
    mp::LinearElasticity,
) where {dim}

    f = zeros(ndofs(dh))
    assembler = start_assemble(K, f)
    nu = getnbasefunctions(cellvalues_u)
    np = getnbasefunctions(cellvalues_p)

    fe = PseudoBlockArray(zeros(nu + np), [nu, np]) # local force vector
    ke = PseudoBlockArray(zeros(nu + np, nu + np), [nu, np], [nu, np]) # local stiffness matrix

    # traction vector
    t = Vec{2}((0.0, -p))
    # cache ɛdev outside the element routine to avoid some unnecessary allocations
    ɛdev = [zero(SymmetricTensor{2,dim}) for i = 1:getnbasefunctions(cellvalues_u)]

    for cell in CellIterator(dh)
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

end;

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

    n_basefuncs_u = getnbasefunctions(cellvalues_u)
    n_basefuncs_p = getnbasefunctions(cellvalues_p)
    u▄, p▄ = 1, 2
    reinit!(cellvalues_u, cell)
    reinit!(cellvalues_p, cell)

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
                Ke[BlockIndex((u▄, u▄), (i, j))] += 2 * mp.G * ɛdev[i] ⊡ ɛdev[j] * dΩ
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
                Ke[BlockIndex((p▄, p▄), (i, j))] += -1 / mp.K * δp * p * dΩ
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
end;

######
#Solve
#######

interpolation_u = Lagrange{2,RefTetrahedron,2}()
interpolation_p = Lagrange{2,RefTetrahedron,1}()


# material
Emod = Eᵣ
Gmod = Emod / 2(1 + ν)
Kmod = Emod * ν / ((1 + ν) * (1 - 2ν))
mp = LinearElasticity(Gmod, Kmod)

# grid, dofhandler, boundary condition
n = 2
grid = create_CMAME1_grid(n, n)
dh = create_dofhandler(grid, interpolation_u, interpolation_p)
dbc = create_bc(dh)

# cellvalues
cellvalues_u, cellvalues_p, facevalues_u = create_values(interpolation_u, interpolation_p)

# assembly and solve
K = create_sparsity_pattern(dh);
K, f = doassemble(cellvalues_u, cellvalues_p, facevalues_u, K, grid, dh, mp);
apply!(K, f, dbc)
u = Symmetric(K) \ f;

# export
filename =
    "CMAME2a" *
    (isa(interpolation_u, Lagrange{2,RefTetrahedron,1}) ? "linear" : "quadratic") *
    "_linear"
vtk_grid(filename, dh) do vtkfile
    vtk_point_data(vtkfile, dh, u)
end

# extract numeric solution
y_points = [Vec((Lᵢₛ, x)) for x in range(0, Lⱼₛ, length=101)];
ph = PointEvalHandler(grid, y_points);
u_points = Ferrite.get_point_values(ph, dh, u, :u);
# extract in x direction
yᵥ = getindex.(u_points, 2)
uyᵥ_num = getindex.(u_points, 2)

uy_anly = 4 * p * Lᵢₛ^3 / (Eᵣ * Lⱼₛ^3)

Plots.plot(
    yᵥ,
    -uyᵥ_num,
    label="Ferrite.jl",
    linecolor=:blue,
    linewidth=4,
    linestyle=:dot,
)
Plots.plot!(
    yᵥ,
    2.4e-3 * ones(length(yᵥ)),
    label="Zerpa 2019 et al",
    linecolor=:green,
    linewidth=4,
    linestyle=:dashdot,
)
Plots.plot!(
    yᵥ,
    uy_anly * ones(length(yᵥ)),
    label="Timoshenko",
    linecolor=:red,
    linewidth=4,
    linestyle=:dashdot,
)
Plots.plot!(
    xlabel=L"y ~ \textrm{[m]}",
    ylabel=L"||u_y(L_i,y)|| ~ \textrm{[m]}",
    legend=:topright,
)
