####################################
#  Plot uniaxial extension results #
####################################
#
# post processing libraries
using Plots, LaTeXStrings
#
# Set to zero values smaller than eps
replace!(x -> x ≤ eps() ? 0 : x, fvalues_apolo_bf)
replace!(x -> x ≤ eps() ? 0 : x, fvalues_closured)
replace!(x -> x ≤ eps() ? 0 : x, fvalues_optim)
if LINEAR_INTENSITY
    replace!(x -> x ≤ eps() ? 0 : x, fvalues_analytic)
end
# load backend
# ------------
gr()
theme(:dracula)
plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2,
    framestyle=:box,
    label=nothing,
    grid=true,
)
# style params
# ------------
colors = (
    analytic=:silver,
    apolo_bf=:green,
    closure=:blue,
    optim=:red,
    gold=:gold,
)
markers = (
    analytic=:hexagon,
    apolo_bf=:none,
    closure=:circle,
    optim=:circle,
    gold=:diamond,
)
lw = 4
linestyles = (
    analytic=:dot,
    apolo_bf=:solid,
    closure=:dash,
    optim=:dash,
    gold=:solid
)
# plot internal apolo solution brute-force
# ----------------------------------------
plot(values(trials_apolo_bf) |> collect, log10.(fvalues_apolo_bf),
    label="Brute-Force Apolo.jl",
    linecolor=colors.apolo_bf,
    linewidth=lw,
    linestyle=linestyles.apolo_bf,
    markershape=markers.apolo_bf,
    markercolor=colors.apolo_bf,
)
# plot Optimizatio.jl solution
# ----------------------------
plot!(collect(values(trials_optim))..., log10.(fvalues_optim),
    label="Optimizations.jl wrapper",
    linecolor=colors.optim,
    linewidth=lw,
    linestyle=linestyles.optim,
    markershape=markers.optim,
    markercolor=colors.optim,
)
# plot closure functional
# ----------------------------
plot!(Evec, log10.(fvalues_closured),
    label="Brute-Force closure function ",
    linecolor=colors.closure,
    linewidth=lw,
    linestyle=linestyles.closure,
    markershape=markers.closure,
    markercolor=colors.closure,
)
if LINEAR_INTENSITY
    # plot analytic functional in the linear intensity case
    # ----------------------------
    plot!(Evec, log10.(fvalues_analytic),
        label="Analytic functional ",
        linecolor=colors.analytic,
        linewidth=lw,
        linestyle=linestyles.analytic,
        markershape=markers.analytic,
        markercolor=colors.analytic,
    )
end
# plot E values
# ----------------------------
vline!([Eᵣ],
    label="gold E = $Eᵣ MPa ",
    linecolor=colors.gold)

vline!([Eᵣ],
    label="Apolo b-f E = $(@sprintf("%.2f", E_value_bf)) MPa",
    linecolor=:transparent)
vline!([Eᵣ],
    label="Optim E = E = $(@sprintf("%.2f", E_value_optim)) MPa",
    linecolor=:transparent)
display(
    plot!(
        xlabel=L"\textrm{Elasticity ~ modulus} ~ E ~ \textrm{[Pa]}",
        ylabel=L"\textrm{Optical-flow ~ functional ~ } F",
        legend=:bottomleft,
    )
)


# Save figure
if LINEAR_INTENSITY
    savefig("./examples/uniaxial_extension/imgs/linear_intensity_results.png")
else
    savefig("./examples/uniaxial_extension/imgs/sinus_intensity_results.png")
end
