####################################
# Brute force resolution algorithm #
####################################

using  InverseProblem:_iterators_unknown_parameters
import InverseProblem: _solve, _initialize!

"""
Brute-force struct algorthim
### Fields:
- `num_params` number of parameters to compute the range
"""
Base.@kwdef struct BruteForceInverseSolver
    num_params::Int = 10


function _solve(invp::MaterialIdentificationProblem, ::BruteForceInverseSolver)

    # compute iterators product
    prams_to_tier = _iterators_unknown_parameters(invp)

    #


end
