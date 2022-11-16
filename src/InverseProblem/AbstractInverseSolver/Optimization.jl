function optimize!(
    ifunctional::AbstractFunctional,
    invp::AbstractInverseProblem;
    alg=BFGS()
)
    # closure over inverse problem
    f = x -> evaluate!(ifunctional, invp, Dict(pᵢ => xᵢ for (pᵢ, xᵢ) in zip(p, x)))
    ofunc = OptimizationFunction(f, Optimization.AutoForwardDiff())

    lb = [search_region[pᵢ][1] for pᵢ in p]
    ub = [search_region[pᵢ][2] for pᵢ in p]
    x0 = (lb + ub) / 2

    prob = OptimizationProblem(ofunc, x0, lb=lb, ub=ub)
    sol = Optimization.solve(prob, alg=alg)
end
