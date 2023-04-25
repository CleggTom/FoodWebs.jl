abstract type AbstractSpecies end


"""
    Species

Type to store niche model parameters for a single species, its thermal optima and a unique ID.  
"""
struct Species <: AbstractSpecies
    n::Float64
    r::Float64
    c::Float64
    Tpk::Float64
    id::UUID
end

"""
    SpeciesParameters

used internaly to move species around ParameterisedCommunites
"""
struct SpeciesParameters
    γ::Float64
    λ::Float64
    μ::Float64
    ϕ::Float64
    ψ::Float64
end

struct ParameterisedSpecies <: AbstractSpecies
    n::Float64
    r::Float64
    c::Float64
    Tpk::Float64
    id::UUID
    p::SpeciesParameters
end


"""
    GeneralisedParameters

Type to contain the parameters for the normalised model.
"""
struct GeneralisedParameters
    #system parameters 
    N::Int64
    n::Vector{Float64}
    A::Matrix{Float64}

    #scale parameters
    α::Vector{Float64}
    β::Matrix{Float64}
    χ::Matrix{Float64}
    ρ::Vector{Float64}
    ρ̃::Vector{Float64}
    σ::Vector{Float64}
    σ̃::Vector{Float64}

    #exponent parameters
    γ::Vector{Float64}
    λ::Matrix{Float64}
    μ::Vector{Float64}
    ϕ::Vector{Float64}
    ψ::Vector{Float64}
end

function get_structural_params(p::GeneralisedParameters)
    return p.β,p.χ,p.ρ,p.ρ̃,p.σ,p.σ̃
end

function get_exponential_params(p::GeneralisedParameters)
    return p.γ, p.λ, p.μ, p.ϕ, p.ψ
end

abstract type  AbstractCommunity end

"""
    Community

Type containing parameters defining a community including its interaction matrix and a vector of species within.
"""
struct Community <: AbstractCommunity
    A::Matrix{Float64}
    sp::Vector{Species}
    ids::Vector{UUID}
    T::Float64
    R::Float64
end

"""
    Community

Type containing parameters defining a community including its interaction matrix and a vector of species within.
"""
struct ParameterisedCommunity <: AbstractCommunity
    A::Matrix{Float64}
    sp::Vector{Species}
    ids::Vector{UUID}
    T::Float64
    R::Float64
    p::GeneralisedParameters
end

"""
    Base.show(io::IO, com::AbstractCommunity)

TBW
"""
Base.show(io::IO, com::AbstractCommunity) = print(io, typeof(com), " N:", length(com.sp)," T:", com.T )


"""
    MetaCommunity

type containing parameters defining a metacommunity including the communties `coms` as well as the matrix `D` which defines "distance" between communties.
"""
struct MetaCommunity
    coms::Array{AbstractCommunity}
    D::Array{Float64}
    T_mat::Array{Float64}
    sp::Vector{Species}
    sp_id::Vector{UUID}
    sp_loc::Dict{UUID,Vector{Int}}
    R::Float64
end

"""
    Base.show(io::IO, mc::MetaCommunity)

TBW
"""
Base.show(io::IO, mc::MetaCommunity) =  print(io, "MetaCommunity M:", length(mc.coms)," ",typeof(mc.coms[1])," Sp: ", length(mc.sp))
