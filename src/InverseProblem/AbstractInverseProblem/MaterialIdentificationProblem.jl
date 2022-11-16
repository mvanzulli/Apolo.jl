
##################################################
# Material Identification Problem implementation #
##################################################

using ..Materials: AbstractParameter
using ..ForwardProblem: AbstractForwardProblem, AbstractForwardProblemSolver
using ..Images: AbstractDataMeasured
using ..InverseProblem: AbstractInverseProblem, AbstractFunctional

export MaterialIdentificationProblem


""" Material identification inverse problem struct.
### Fields:
- `fproblem`    -- forward problem formulation.
- `fsolver`     -- forward problem solver.
- `datam`       -- measured data input to solve the inverse problem.
- `f`           -- target function or functional to be minimized.
- `roi`         -- region of interest (in space) where the functional is going to be evaluated.
- `sregion`     -- region where the parameters are going to be explored.
"""
struct MaterialIdentificationProblem{
    FP<:AbstractForwardProblem,
    FSOL<:AbstractForwardProblemSolver,
    DM<:AbstractDataMeasured,
    F<:AbstractFunctional,
    R,
    P<:AbstractParameter,
    PRANGE<:AbstractVector} <: AbstractInverseProblem
    #Fields:
    fproblem::FP
    fsolver::FSOL
    datam::DM
    f::F
    roi::R
    sregion::Dict{P,PRANGE}

end

"Constructor without a search region defined will consider the feasible region."
function MaterialIdentificationProblem(
    fproblem::AbstractForwardProblem,
    fsolver::AbstractForwardProblemSolver,
    datam::AbstractDataMeasured,
    func::AbstractFunctional,
    roi::Function,
)

    uparams = unknown_parameters(fproblem)
    search_default_region = Dict{AbstractParameter, AbstractVector}()
    for uparam in uparams
        search_default_region[uparam] = range(uparam)
    end

    return MaterialIdentificationProblem(
        fproblem, fsolver, datam, func, roi, search_default_region
        )

end
