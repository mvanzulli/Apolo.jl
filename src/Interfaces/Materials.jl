#######################################################
# Main types and function sto solve a Forward Problem  #
#######################################################

"""
Module defining the materials interface. Each material consists of a data type with one or more `Parameter` fields.
"""
module Materials

using AutoHashEquals: @auto_hash_equals


export AbstractMaterial, ConstitutiveParameter, SVK
export material, setmaterial!, has_material, label, value, setval!, feasible_region,
    has_feasible_region, set_feasible_region!, parameter_index, set!, setval!, model, parameters

const REALS = (-Inf, Inf) # Tuple to check if (v₁, v₂) ∈ (-∞, +∞)
const DEFAULT_NUMBER_OF_PARAMS_RANGE = 10 # Default number of parameters when a range is created in a LinRange

""" Abstract supertype for material parameters.

The following methods are provided by the interface:

- `material(p)`                     -- returns the parmeter's material.
- `has_material(p)`                 -- returns `true` if material of the parameter is assigned.
- `label(p)`                        -- returns the label of the parameter `p`.
- `ismissing(p)`                    -- returns `true` if the value of parameter `p` is missing.
- `value(p)`                        -- returns the parameter `p` value.
- `range(p)`                        -- returns the parameter `p` range.
- `has_feasible_region(p)`          -- returns `true` if the parameter `p` has a constrained range defined.
- `feasible_region(p)`              -- returns the feasible region for the parameter `p`.
- `setval!(p, val)`                 -- sets the value `val` to the parameter `p`.
- `set_feasible_region!(p, range)`  -- sets the limits (`p₁`, `p₂`) to the parameter `p`.
- `setmaterial!(p, mlabel)`         -- sets a material label for a the paretmer.
"""
abstract type AbstractParameter end

"Checks if the parameter `p` has a constrainted region defined."
has_feasible_region(p::AbstractParameter) = any(!, isinf.(feasible_region(p)))

"Returns the paramter `p`` bounds."
feasible_region(p::AbstractParameter) = p.fregion

"Returns the label of the parameter `p`."
label(p::AbstractParameter) = Symbol(p.label)

"Returns the material label of the parameter `p`."
material(p::AbstractParameter) = Symbol(p.material)

"Checks if the parameter `p` has a constrainted range defined."
has_material(p::AbstractParameter) = material(p) ≠ :no_assigned

"Returns `true` if label of the parameter `p` and `plabel` match."
_islabel(p::AbstractParameter, plabel::Symbol) = label(p) == plabel

"Returns `true` if the value of parameter `p` is missing."
Base.ismissing(p::AbstractParameter) = ismissing(p.val)

"Returns the parameter `p` value."
value(p::AbstractParameter) =
    !ismissing(p) ? p.val : throw(ArgumentError("The value of $(label(p)) is unknown"))

"Return the parameter `p` range."
Base.range(p::AbstractParameter,
    num_p::Int=DEFAULT_NUMBER_OF_PARAMS_RANGE) =
    !has_feasible_region(p) ?
    throw(
        ArgumentError(
            "The range is not specified for $(label(p))"
        )
    ) : LinRange(feasible_region(p)[1], feasible_region(p)[2], num_p)

"Checks if a value is into the range of feasible parameters."
Base.:∈(val::Number, p::AbstractParameter) = feasible_region(p)[1] ≤ val ≤ feasible_region(p)[2]

"Checks if a value is outside the range of feasible parameters."
Base.:∉(val::Number, p::AbstractParameter) = !(val ∈ p)

"Sets the value `val` to the parameter `p`."
function setval!(p::AbstractParameter, val::Number)

    val ∉ p && throw(ArgumentError(
        "The value $val is not inside the parameter range = $(feasible_region(p))"
        ))

    p.val = val

    return p

end

"Sets the range (`pₘᵢₙ`, `pₘₐₓ`) to the parameter `p`."
function set_feasible_region!(p::AbstractParameter, pₘᵢₙ::Number, pₘₐₓ::Number)

    has_feasible_region(p) && @warn(
        "The feasible region of p = $p is being modified ($(feasible_region(p)) => ($pₘᵢₙ, $pₘₐₓ))"
        )
    p.fregion = (pₘᵢₙ, pₘₐₓ)

    return p

end

"Sets the material label `mlabel` to the parameter `p`."
function setmaterial!(p::AbstractParameter, mlabel::Symbol)

    has_material(p) && @warn("The paramter's material has been overwritten")
    p.material = mlabel

    return p

end

""" Constitutive parameter struct.
### Fields:

- `label`             -- label
- `val`               -- value
- `feasible_region`   -- parameter guess feasible_region
- `material`          -- material symbol

"""
@auto_hash_equals mutable struct ConstitutiveParameter <: AbstractParameter
    label::Symbol
    val::Union{Number,Missing}
    fregion::Union{NTuple{2,Number},Missing}
    material::Symbol
    function ConstitutiveParameter(label, val=missing, fregion=REALS, mat=:no_assigned)
        new(Symbol(label), val, fregion, Symbol(mat))
    end
end

""" Abstract supertype for a material.

The following methods are provided by the interface:

- `model(m)`              -- return a string with the material model (defaults to the materials' type label)
- `parameters(m)`         -- return a tuple with the material parameters (defaults to the materials' fields which are of type `Parameter`)
- `label(m)`           -- return material label
- `getindex(m, plabel)`    -- return the index of the parameter with label `plabel` into the list of material `m` parameters.
- `get(m, p)` or `m[:p]`  -- return the parameter `p` of the material `m`.
- `value(m, plabel)`      -- return the value of the material parameter with label `plabel`.
- `range(m, plabel)`    -- return the range of the material parameter with label `plabel`.
- `set!(m, p)`            -- set parameter `p` into material `m`.
- `setval!(m, p, pval)`   -- set value `pval` to parameter labeld `plabel` into material `m`.
- `range(p)`           -- return the parameter `p` range.
- `set_val(p, val)`       -- set the val `val` to the parameter `p`.
- `set_range(p, range)`   -- set the range `range` to the parameter `p`.

"""
abstract type AbstractMaterial end

model(::Type{T}) where {T<:AbstractMaterial} = string(T)

function parameters(m::T) where {T<:AbstractMaterial}
    Tuple([getfield(f, n) for n in fieldlabels(T) if fieldtype(T, n) isa Parameter])
end

"Returns label."
label(m::AbstractMaterial) = ""

"Returns the index of the parameter with label `plabel` into the list of material `m` parameters."
function parameter_index(m::AbstractMaterial, plabel::Symbol)
    index_param = findall(_islabel.(parameters(m), plabel))
    length(index_param) !== 1 ?
    throw(ArgumentError("parameter $plabel ∉ mat $m or is not unique")) : first(index_param)
end

"Returns the parameter `p` of the material `m`."
Base.getindex(m::AbstractMaterial, plabel::Symbol) = getfield(m, plabel)

"Returns the value of the material parameter with label `plabel`."
value(m::AbstractMaterial, plabel::Symbol) = value(getindex(m, plabel))

"Returns the range of the material parameter with label `plabel`."
Base.range(m::AbstractMaterial, plabel::Symbol) = range(m[plabel])

" Set parameter `p` into material `m`."
function Base.replace!(m::AbstractMaterial, pair::Pair{Symbol,<:AbstractParameter})

    @assert ismutable(m) && hasfield(typeof(m), pair[1])
    mat_name = Symbol(label(m))
    mat_name != material(pair[2]) && setmaterial!(pair[2], mat_name)
    setfield!(m, pair[1], pair[2])

    return m

end

" Set value `pval` to parameter labeld `plabel` into a material `m`."
function setval!(m::AbstractMaterial, plabel::Symbol, pval::Number)

    @assert hasfield(typeof(m), plabel)
    setval!(m[plabel], pval)

    return m
end

"Returns [`true`](@ref) if the parameter with label `plabel` is [`missing`](@ref) into material `m`."
Base.ismissing(m::AbstractMaterial, plabel::Symbol) = ismissing(m[plabel])

# ==============================
# Concrete implementations
# ==============================

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
"Constructor with a string label for SVK amaterial."
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


end # end module
