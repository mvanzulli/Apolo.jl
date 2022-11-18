# Plot min max and gold solutions
# ------------------------------
# general plot params
const LINE_WIDTH = 4
const MARKER_SIZE = 3
# min,max,gold
const GOLD_COLOR = :goldenrod3
const GOLD_MARKER = :circle
const MIN_COLOR = :darkorange3
const MIN_MARKER = :hexagon
const MAX_COLOR = :silver
const MAX_MARKER = :rect

" Plots max, min and gold solutions"
function plot_extremaΩ_sols(
    Eᵥ::Vector,νᵥ::νT,
    gold_solution::ForwardProblemSolution,
    min_solution::ForwardProblemSolution,
    max_solution::ForwardProblemSolution,
    num_points_border = 20,
    ) where νT<:Union{Vector,LinRange}

        # Check Eᵥ contains only one value
        @assert length(Eᵥ) == 1 "In this method Eᵥ must contain only one value"

        # plot min solution|
        Ωₘᵢₙ = SolutionBoundary(min_solution, num_points_border)
        scatter(Tuple.(get_elements(Ωₘᵢₙ)),
                label = "min sol (E,ν) = ($(Eᵥ[1]),$(νᵥ[1]))",
                color = MIN_COLOR,
                markershape = MIN_MARKER,
                markersize = MARKER_SIZE,
                )

        # plot max solution
        Ωₘₐₓ = SolutionBoundary(max_solution, num_points_border)
        scatter!(Tuple.(get_elements(Ωₘₐₓ)),
            label = "max sol (E,ν) = ($(Eᵥ[end]),$(νᵥ[end]))",
            color = MAX_COLOR,
            markershape = MAX_MARKER,
            markersize = MARKER_SIZE)

        # plot gold solution
        Ωgold = SolutionBoundary(gold_solution, num_points_border)
        scatter!(Tuple.(get_elements(Ωgold)),
        label = "gold sol (E,ν) = ($(Eᵣ),$(νᵣ))",
        color = GOLD_COLOR,
        markershape = GOLD_MARKER,
        markersize = MARKER_SIZE,)

        # display plot
        display(scatter!(
            xlabel=L"x ~ \textrm{[m]}",
            ylabel=L"y ~ \textrm{[m]}",
            legend=:best)
            )

end
#
# plot them
plot_extremaΩ_sols(Eᵥ, νᵥ, gold_solution, min_solution, max_solution)
# Plot jaccard number
# ------------------------------
# plot params
const JACCARD_COLOR = :deepskyblue3
const JACCARD_MARKER = :circle
"Plots jaccards number when ν is changing"
function plot_jaccard_ν(νᵥ::νT, jaccards) where νT<:Union{Vector,LinRange}
    # Plot Jaccard's number
    plot(
        νᵥ,
        jaccards[1,:],
        label="Jaccard's number",
        linecolor=JACCARD_COLOR,
        markercolor=JACCARD_COLOR,
        linewidth=LINE_WIDTH,
        linestyle=:solid,
        markershape=JACCARD_MARKER,
        markersize=MARKER_SIZE,
    )
   # Plot gold ν
    vline!(
        [νᵣ],
        label="νᵣ = $νᵣ",
        linecolor=GOLD_COLOR,
        linewidth=LINE_WIDTH,
        linestyle=:dot,
        markershape=GOLD_MARKER,
        markercolor=GOLD_COLOR,
        markersize=MARKER_SIZE,
    )
    # add legends and labels
    display(
        plot!(
            xlabel=L"\nu",
            ylabel=L"J",
            legend=:best,
        )
    )
end
#
plot_jaccard_ν(νᵥ , jaccards)
# Plot Hausdorff distance
# ------------------------------
# plot params
const HAUSDORFF_COLOR = :purple2
const HAUSDORFF_MARKER = :circle
function plot_hausdorff_ν(νᵥ::νT, hausdorffs) where νT<:Union{Vector,LinRange}
    # Plot Jaccard's number
    plot(
        νᵥ,
        hausdorffs[1,:],
        label="Hausdorff's number",
        linecolor=HAUSDORFF_COLOR,
        linewidth=LINE_WIDTH,
        linestyle=:solid,
        markershape=HAUSDORFF_MARKER,
        markercolor=HAUSDORFF_COLOR,
        markersize=MARKER_SIZE,
    )
   # Plot gold ν
    vline!(
        [νᵣ],
        label="νᵣ = $νᵣ",
        linecolor=GOLD_COLOR,
        linewidth=LINE_WIDTH,
        linestyle=:dot,
        markershape=GOLD_MARKER,
        markercolor=GOLD_COLOR,
        markersize=MARKER_SIZE,
    )
    # add legends and labels
    display(
        plot!(
            xlabel=L"\nu",
            ylabel=L"H",
            legend=:best,
        )
    )
end
#
plot_hausdorff_ν(νᵥ, hausdorffs)
# Plot geometric similarty function
# -----------------------------------
# plot params
const GSF_COLOR = :seagreen
const HAUSDORFF_MARKER = :circle
function plot_gsf_ν(νᵥ::νT, gsfs) where νT<:Union{Vector,LinRange}
    # Plot Jaccard's number
    plot(
        νᵥ,
        gsfs[1,:],
        label="GSF number",
        linecolor=GSF_COLOR,
        linewidth=LINE_WIDTH,
        linestyle=:solid,
        markershape=HAUSDORFF_MARKER,
        markercolor=GSF_COLOR,
        markersize=MARKER_SIZE,
    )
   # Plot gold ν
    vline!(
        [νᵣ],
        label="νᵣ = $νᵣ",
        linecolor=GOLD_COLOR,
        linewidth=LINE_WIDTH,
        linestyle=:dot,
        markershape=GOLD_MARKER,
        markercolor=GOLD_COLOR,
        markersize=MARKER_SIZE,
    )
    # add legends and labels
    display(
        plot!(
            xlabel=L"\nu",
            ylabel=L"GSF",
            legend=:best,
        )
    )
end
#
plot_gsf_ν(νᵥ, gsfs)
