"""
Module defining linear elastic materials.
"""
module LinearElastic

using AutoHashEquals: @auto_hash_equals

using ..Materials: AbstractMaterial, ConstitutiveParameter
using ..Materials: setmaterial!, value

import ..Materials: label, parameters, model

export SVK, lamé_params

""" SVK material struct.

### Fields:

`E` -- Elasticity modulus.
`ν` -- Poisson's ratio.
`label` -- Label to recognize material

"""
@auto_hash_equals mutable struct SVK <: AbstractMaterial
    E::ConstitutiveParameter
    ν::ConstitutiveParameter
    label::Symbol
    function SVK(E, ν, label=:no_assigned)

        setmaterial!(E, label)
        setmaterial!(ν, label)

        return new(E, ν, label)
    end
end

"Constructor with a string label for SVK material."
function SVK(E, ν, label::String)

    label = Symbol(label)
    setmaterial!(E, label)
    setmaterial!(ν, label)

    return SVK(E, ν, label)
end

"Returns SVK material model label."
model(::SVK) = "SVK"

"Returns SVK parameters tuple."
parameters(m::SVK) = (m.E, m.ν)

"Returns SVK label"
label(m::SVK) = string(m.label)

" Extract svk material parameters to use with ferrite nomenclature"
function lamé_params(svk::SVK)

    E = value(svk[:E])
    ν = value(svk[:ν])

    # Compute Lamé parameters λ and μ (μ = G)
    μ = E / 2(1 + ν)
    λ = E * ν / ((1 + ν) * (1 - 2ν))

    return μ, λ
end

end
