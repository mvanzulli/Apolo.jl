####################################
# Brute force resolution algorithm #
####################################

using ..InverseProblem: InverseProblemSolution
using ..InverseProblem: functional, values, trials, evaluate!, _set_optim_done!
using ..Utils: ScalarWrapper
import ..InverseProblem: solve

export BruteForce, nparams_foreach_param

"Returns the unknown parameters iterators of a given inverse problem `iproblem`."
function _iterators_unknown_parameters(iproblem::AbstractInverseProblem)

    uparams = unknown_parameters(iproblem)
    sregion = search_region(iproblem)
    sregion_uparams = [sregion[u] for u in uparams]

    # combine them into a single vector
    iters = vec([p_combination for p_combination in Iterators.product(sregion_uparams...) ])
    # set parameters dicts iterators
    set_params_iters = Vector{Dict}()
    for iter_p in iters
        d = Dict{AbstractParameter,eltype(eltype(iters))}()
        for (ip, p) in enumerate(uparams)
            d[p] = iter_p[ip]
        end
        push!(set_params_iters, d)
    end

    return set_params_iters
end

"Returns the unknown parameters iterators of a given inverse problem `iproblem` and a number of parameters."
function _iterators_unknown_parameters(iproblem::AbstractInverseProblem, num_params::Int)

    uparams = unknown_parameters(iproblem)
    sregion = search_region(iproblem)
    params_ranges_default = [sregion[u] for u in uparams]

    # Update the region with same limits and different length
    sregion_uparams = [
        LinRange(params_ranges_default[i].start, params_ranges_default[i].stop, num_params)
        for i in eachindex(params_ranges_default)]

    # combine them into a single vector
    iters = vec([p_combination for p_combination in Iterators.product(sregion_uparams...) ])

    # set parameters dicts iterators
    set_params_iters = Vector{Dict}()
    for iter_p in iters
        d = Dict{AbstractParameter,eltype(eltype(iters))}()
        for (ip, p) in enumerate(uparams)
            d[p] = iter_p[ip]
        end
        push!(set_params_iters, d)
    end

    return set_params_iters

end


"""
Brute-force struct algorthim
### Fields:
- `num_params` -- number of paramaters that are tried for each unknown parameter.
"""
Base.@kwdef struct BruteForce <:AbstractInverseProblemSolver
    nparam::Int = 10
    optim_done::ScalarWrapper = ScalarWrapper(false)
end

"Brute force algorithm constructor with a `nparmas` number of parameters "
BruteForce(naparams::Int) = BruteForce(naparams, ScalarWrapper(false))

"Solves the material identification problem."
function solve(invp::MaterialIdentificationProblem, bfs::BruteForce)

    #compute parameters set of iterators
    prams_to_tier = _iterators_unknown_parameters(invp, nparams_foreach_param(bfs))

    # exhaustive evaluation
    func = functional(invp)
    [evaluate!(func, invp, params) for params in prams_to_tier]


    _set_optim_done!(bfs)

    #
    # Main.@infiltrate
    return InverseProblemSolution(invp, func, bfs)
end

"Returns the number of paramaters that are tried for each unknown parameter."
nparams_foreach_param(bfs::BruteForce) = bfs.nparam
