#######################################################
# Main types and function sto solve a Forward Problem  #
#######################################################

#TODO: add Iterator for materials for params in p

"""
Module defining the materials interface. Each material consists of a data type with one or more `Parameter` fields.
"""
module Materials

# Add libraries to use
using AutoHashEquals: @auto_hash_equals

# Export interface functions and types
export AbstractMaterial, Parameter, SVK

export name, isname, label, value, setval!, hasrange,
    setrange!, parameter_index, set!, setval!, model, parameters

""" Abstract supertype for material parameters.

The following methods are provided by the interface:

- `name(p)`                 -- return the name of the parameter `p`.
- `isname(p, pname)`        -- return `true` if name of the parameter `p` and `pname` match.
- `ismissing(p)`            -- return `true` if the value of parameter `p` is missing.
- `value(p)`               -- return the parameter `p` value.
- `range(p)`             -- return the parameter `p` range.
- `setval!(p, val)`         -- set the value `val` to the parameter `p`.
- `setrange!(p, range)`     -- set the range (`p₁`, `p₂`) to the parameter `p`
"""
abstract type AbstractParameter end

# Constant variables or types
const REALS = (-Inf, Inf) # Tuple to check if (v₁, v₂) ∈ (-∞, +∞)
const DEFAULT_NUMBER_OF_PARAMS_RANGE = 10 # Default number of parameters when a range is created in a LinRange
""" Parameter struct.
### Fields:

- `name`     -- name
- `val`      -- value
- `range`    -- parameter guess range

"""
@auto_hash_equals mutable struct Parameter <: AbstractParameter
    name::Symbol
    val::Union{Number,Missing}
    range::Union{NTuple{2,Number},Missing}
    function Parameter(name, val=missing, range=REALS)
        new(Symbol(name), val, range)
    end
end

## Methods for abstract parameters:
" Return the name of the parameter `p`. "
name(p::AbstractParameter) = Symbol(p.name)

" Return `true` if name of the parameter `p` and `pname` match. "
isname(p::AbstractParameter, pname::Symbol) = name(p) == pname

" Return `true` if the value of parameter `p` is missing. "
Base.ismissing(p::AbstractParameter) = ismissing(p.val)

" Return the parameter `p` value. "
value(p::AbstractParameter) =
    !ismissing(p) ? p.val : throw(ArgumentError("The value of $(name(p)) is unknown"))

" Checks if the parameter `p` has a constrainted range defined. "
hasrange(p::AbstractParameter) = any(!, isinf.(p.range))

"Return the parameter `p` range."
Base.range(p::AbstractParameter,
    num_p::Int=DEFAULT_NUMBER_OF_PARAMS_RANGE) =
    !hasrange(p) ?
    throw(
        ArgumentError(
            "The range is not specified for $(name(p))"
        )
    ) : LinRange(p.range[1], p.range[2], num_p)

"Checks if a value is into the range of feaseable parameters."
Base.:∈(val::Number, p::Parameter) = p.range[1] ≤ val ≤ p.range[2]

"Checks if a value is outside the range of feaseable parameters."
Base.:∉(val::Number, p::Parameter) = !(val ∈ p)

"Sets the value `val` to the parameter `p`."
function setval!(p::AbstractParameter, val::Number)
    val ∉ p && throw(ArgumentError("The value $val is not inside the parameter range = $(p.range)"))
    p.val = val
    return p
end

"Sets the range (`pₘᵢₙ`, `pₘₐₓ`) to the parameter `p`."
function setrange!(p::AbstractParameter, pₘᵢₙ::Number, pₘₐₓ::Number)
    p.range = (pₘᵢₙ, pₘₐₓ)
    return p
end

""" Abstract supertype for a material.

The following methods are provided by the interface:

- `model(m)`              -- return a string with the material model (defaults to the materials' type name)
- `parameters(m)`         -- return a tuple with the material parameters (defaults to the materials' fields which are of type `Parameter`)
- `label(m)`           -- return material label
- `getindex(m, pname)`    -- return the index of the parameter with name `pname` into the list of material `m` parameters.
- `get(m, p)` or `m[:p]`  -- return the parameter `p` of the material `m`.
- `value(m, pname)`      -- return the value of the material parameter with name `pname`.
- `range(m, pname)`    -- return the range of the material parameter with name `pname`.
- `set!(m, p)`            -- set parameter `p` into material `m`.
- `setval!(m, p, pval)`   -- set value `pval` to parameter named `pname` into material `m`.
- `range(p)`           -- return the parameter `p` range.
- `set_val(p, val)`       -- set the val `val` to the parameter `p`.
- `set_range(p, range)`   -- set the range `range` to the parameter `p`.

"""
abstract type AbstractMaterial end

model(::Type{T}) where {T<:AbstractMaterial} = string(T)

function parameters(m::T) where {T<:AbstractMaterial}
    Tuple([getfield(f, n) for n in fieldnames(T) if fieldtype(T, n) isa Parameter])
end

" Return label. "
label(m::AbstractMaterial) = ""

" Return the index of the parameter with name `pname` into the list of material `m` parameters. "
function parameter_index(m::AbstractMaterial, pname::Symbol)
    index_param = findall(isname.(parameters(m), pname))
    length(index_param) !== 1 ?
    throw(ArgumentError("parameter $pname ∉ mat $m or is not unique")) : first(index_param)
end

" Return the parameter `p` of the material `m`. "
Base.getindex(m::AbstractMaterial, pname::Symbol) = getfield(m, pname)

" Return the value of the material parameter with name `pname`. "
value(m::AbstractMaterial, pname::Symbol) = value(getindex(m, pname))

" Return the range of the material parameter with name `pname`. "
Base.range(m::AbstractMaterial, pname::Symbol) = range(m[pname])

" Set parameter `p` into material `m`. "
function Base.replace!(m::AbstractMaterial, pair::Pair{Symbol,<:AbstractParameter})
    @assert ismutable(m) && hasfield(typeof(m), pair[1])
    setfield!(m, pair[1], pair[2])
    return m
end

" Set value `pval` to parameter named `pname` into a material `m`. "
function setval!(m::AbstractMaterial, pname::Symbol, pval::Number)
    @assert hasfield(typeof(m), pname)
    setval!(m[pname], pval)
    return m
end

" Return [`true`](@ref) if the parameter with name `pname` is [`missing`](@ref) into material `m`. "
Base.ismissing(m::AbstractMaterial, pname::Symbol) = ismissing(m[pname])

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
    E::Parameter
    ν::Parameter
    label::Symbol
    function SVK(E, ν, label="")
        return new(E, ν, Symbol(label))
    end
end
"Returns SVK material model name."
model(::SVK) = "SVK"

"Returns SVK parameters tuple."
parameters(m::SVK) = (m.E, m.ν)

"Returns SVK label"
label(m::SVK) = string(m.label)

" Extract svk material parameters to use with ferrite nomenclature"
function lamé_params(svk::SVK)
    # Extract E and ν parameters
    E = value(svk[:E])
    ν = value(svk[:ν])

    # Compute Lamé parameters λ and μ (μ = G)
    μ = E / 2(1 + ν)
    λ = E * ν / ((1 + ν) * (1 - 2ν))

    return μ, λ
end


end # end module
