function structural_parameters(A)
    N = size(A)[1]
    
    β = zeros(N,N)
    χ = zeros(N,N)
    for i = 1:N
        for j = 1:N
            #relative loss of j due to predation by i
            β[i,j] = A[i,j] / sum(A[:,j]) 
            #relative contributuion of j to total predation by i
            χ[i,j] = A[i,j] / sum(A[i,:]) 
        end
    end

    β[isnan.(β)] .= 0.0
    χ[isnan.(χ)] .= 0.0

    #total growth from produciton
    ρ = 1 .* (sum(A,dims=2) .!= 0)[:]
    #growth from predation
    ρ̃ = 1 .- ρ
    
    #loss from predation
    σ = 1 .* (sum(A,dims=1) .!= 0)[:]
    #loss from mortality
    σ̃ = 1 .- σ

    return(β,χ,ρ,ρ̃,σ,σ̃)
end

function random_parameters(N)
    #exponent
    γ = rand(N) .+ 0.5 #[0.5,1.5]
    λ = ones(N,N)
    μ = rand(N) .+ 1.0 #[1.0,2.0]
    ϕ = rand(N)
    ψ = rand(N) .+ 0.5

    return(γ,λ,μ,ϕ,ψ)
end

"""
    generalised_parameters(A,n)

Randomly samples a set of generalised parameters from a given interaction matrix and set of niches 
"""
function generalised_parameters(A,n; R::Float64 = 42.0)
    N = size(A)[1]
    #scale
    α = (R .^ n) .^ (-0.25)

    β,χ,ρ,ρ̃,σ,σ̃ = structural_parameters(A)
    γ,λ,μ,ϕ,ψ = random_parameters(N)
    

    return GeneralisedParameters(N,n,A,α,β,χ,ρ,ρ̃,σ,σ̃,γ,λ,μ,ϕ,ψ)
end


"""
    generalised_parameters(com::Community)

generate parameters directly from a Community object
"""
function generalised_parameters(com::Community)
    generalised_parameters(com.A, [s.n for s = com.sp], R = com.R)
end


"""
    generalised_jacobian!(J,p::GeneralisedParameters)

Constructs the jaccobian in place as defined by the generalised model. Takes parameter values (including network structure) from the `GeneralisedParameters` object `p`.
"""
function generalised_jacobian!(J,p::GeneralisedParameters)
    @assert size(J) == size(p.A)
    #construct generalised_jacobian
    for i = 1:p.N
        for j = 1:p.N
            if i == j #intra
                J[i,j] = p.ρ̃[i] * p.ϕ[i] + #production
                         p.ρ[i] * (p.γ[i] * p.χ[i,j] * p.λ[i,j] + p.ψ[i]) - #predation
                         p.σ̃[i] * p.μ[i] #mortality
                for k = 1:p.N
                    J[i,j] -= p.σ[i] * p.β[k,i] * p.λ[k,i] * ((p.γ[k] - 1) * p.χ[k,i] + 1) #predation
                end
            else #inter
                J[i,j] = p.ρ[i] * p.γ[i] * p.χ[i,j] * p.λ[i,j] - #predation
                         p.σ[i] * p.β[j,i] * p.ψ[j]
                for k = 1:p.N
                    J[i,j] -= p.σ[i] * (p.β[k,i] * p.λ[k,j] * (p.γ[k] - 1) * p.χ[k,j])
                end
            end
            J[i,j] *= p.α[i]
        end
    end
    
end

"""
    generalised_jacobian(p::GeneralisedParameters)

generate jacobian from parameter set p
"""
function generalised_jacobian(p::GeneralisedParameters)
    J = zeros(p.N, p.N)
    generalised_jacobian!(J,p)
    return J 
end


#Calculate stablilty
"""
    max_real_eigval(J)

Get the real part of the maximum eigenvalue from J allowing for complex values. 
"""
function max_real_eigval(J)
    λ = eigvals(J)

    if length(λ)>0
        if isa(λ[end], Complex)
            return(λ[end].re)
        else
            return(λ[end])
        end
    end

    return(NaN)
end


"""
    communtiy_stability(J,p::GeneralisedParameters)

Get the real part of the leading eigenvalue from parameter set p. J is provided to avoid allocation
"""
function communtiy_stability!(J,p::GeneralisedParameters)
    if p.N == 0 
        return NaN
    end

    generalised_jacobian!(J,p::GeneralisedParameters)
    return max_real_eigval(J) 
end

"""
    communtiy_stability(p::GeneralisedParameters)

Non-inplace version.
"""
function communtiy_stability(p::GeneralisedParameters)
    J = zeros(size(p.A))
    return communtiy_stability!(J,p::GeneralisedParameters)
end

"""
    communtiy_stability(J, c::Community)

Calculates real part of the leading eigenvalue directly from a Communtiy object. 
"""
function communtiy_stability!(J, c::Community)
    #make params
    p = generalised_parameters(c)
    return communtiy_stability!(J, p)
end

"""
    communtiy_stability(c::Community)

Calculates real part of the leading eigenvalue directly from a Communitiy object. 
"""
function communtiy_stability(c::Community)
    #make params
    J = zeros(size(c.A))
    return communtiy_stability!(J, c)
end

function communtiy_stability!(J, c::ParameterisedCommunity)
    communtiy_stability!(J, c.p)
end

function communtiy_stability(c::ParameterisedCommunity)
    communtiy_stability(c.p)
end


function proportion_stable_webs(J,c::Community; N_trials::Int = 100)
    sum([communtiy_stability(J, c) for i = 1:N_trials] .< 0) / N_trials
end

"""
    proportion_stable_webs(c::Community, N::Int = 100)

Calculates the proportion stable parameter configurations for a Community.   
"""
function proportion_stable_webs(c::Community; N_trials::Int = 100)
    J = zeros(size(c.A))
    sum([communtiy_stability(J, c) for i = 1:N_trials] .< 0) / N_trials
end
 
"""
    proportion_stable_webs(mc::MetaCommunity; N_trials::Int = 100)

Calculates proprtion stable webs across a metacommuntiy, returning an array of proportions. 
"""
function proportion_stable_webs(mc::MetaCommunity; N_trials::Int = 100)
    props = zeros(size(mc.coms))
    for (i,c) = enumerate(mc.coms)
        if length(c.sp) > 0
            J = zeros(size(c.A))
            props[i] = proportion_stable_webs(J,c, N_trials = N_trials)
        end
    end
    return props 
end


function stable_parameterisation(com::Community, N_trials = 100)
    k = 1
    J = zeros(size(com.A))
    while k < N_trials
 
        p = generalised_parameters(com)
        if communtiy_stability!(J, p) < 0
            return parameterised_community(com, p)
        end
        k += 1
    end
    return com
end

function time2stable(com)
    t_start = time()

    if stable_parameterisation(com, Inf) == false
        return Inf
    end

    t_end = time()

    return t_end - t_start
end

function metacommuntiy_stability(mc)
    [communtiy_stability(c) for c = mc.coms]
end
