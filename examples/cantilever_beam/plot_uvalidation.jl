# -------------------------------------------------------------------------------------------------
# This is the Ferrite script of a one dimensional identification problem of Zerpa 2019 usin IdenGPU
# ------------------------------------------------------------------------------------------------
# u vals
u = sol.valdofs
# extract numeric solution
# create a point handler
ph = PointEvalHandler(grid, y_points)
# eval
u_points = Ferrite.get_point_values(ph, dh, u, :u);
# displacement in y direction
yᵥ = getindex.(u_points, 2)
uyᵥ_num = getindex.(u_points, 2)
# analytic solution
uy_anly = 4 * pₓ * Lᵢₛ^3 / (Eᵣ * Lⱼₛ^3)
# Plot uy at the righ
# ------------------------------
gegetfaceset(grid, "top")

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
    label="Euler-Bernoulli beam",
    linecolor=:red,
    linewidth=4,
    linestyle=:dashdot,
)
Plots.plot!(
    xlabel=L"y ~ \textrm{[m]}",
    ylabel=L"||u_y(L_i,y)|| ~ \textrm{[m]}",
    legend=:topright,
)