
"""
    species(C::Float64)

Constructor function for Sp type. Randomly generates niche parameters given a connectance value
"""
function species(C::Float64, Tpk::Float64; n::Float64 = rand())
    #niche params
    β = (1 / (2*C)) - 1
    # n = rand()
    r = n * rand(Beta(1.0, β))
    c = rand(Uniform(r / 2, n))
    #extra params
    uuid = UUIDs.uuid1()
    return Species(n,r,c,Tpk,uuid)
end

species(C::Float64; n::Float64 = rand()) = species(C, rand(), n = n)

"""
    species(C::Float64, Tpk::Vector{Float64}, N::Int64)

efficent method for generating multiple species at a time.
"""
function species(C::Float64, Tpk::Vector{Float64}, N::Int64; n::Vector{Float64} = rand(N))
    β = (1 / (2*C)) - 1
    # n = rand(N)
    r = n .* rand(Beta(1.0, β), N)
    c = rand.(Uniform.(r / 2, n))
    uuid = [UUIDs.uuid1() for i = 1:N]

    return Species.(n,r,c,Tpk,uuid)
end

species(C::Float64, N::Int64; n::Vector{Float64} = rand(N)) = species(C, rand(N), N, n = n)


"""
    parameterised_species(sp::Species, p::SpeciesParameters)

create parameterised speceies with given expoential parameters
"""
function parameterised_species(sp::Species, p::SpeciesParameters)
    ParameterisedSpecies(sp.n, sp.r, sp.c, sp.Tpk, sp.id, p)
end

"""
    parameterised_species(com::Community, id::UUID)

extract a parameterised species from a ParameterisedCommunity
"""
function parameterised_species(com::ParameterisedCommunity, id::UUID)
    #make param set
    indx = findall(id .== com.ids)

    if length(indx) > 0
        if length(indx) == 1
            indx = indx[1]
            sp_p = SpeciesParameters(com.p.γ[indx],com.p.λ[indx,indx],com.p.μ[indx],com.p.ϕ[indx], com.p.ψ[indx])
            return parameterised_species(com.sp[indx], sp_p)
        else
            error("mulitple sp with id in the community")
        end
    else
        error("species with id not in community")
    end
    
end

"""
    isequalcommunity(a::AbstractCommunity, b::AbstractCommunity)

Test if two community objects are equal. This forces comparison by equality (`==`) as opposed to identity ('===') which the default `isequal` method uses. 
"""
function isequalcommunity(a::Community, b::Community)
    typeof(a) != typeof(b) && return(false)
    #get keys
    com_keys = fieldnames(Community)
    #test if all fields are the same
    return getfield.(Ref(a), com_keys) == getfield.(Ref(b), com_keys) 
end

function isequalcommunity(a::ParameterisedCommunity, b::ParameterisedCommunity)
    typeof(a) != typeof(b) && return(false)
    #get keys
    #test if all fields are the same
    for k = fieldnames(ParameterisedCommunity)
        if k == :p
            for kp = fieldnames(GeneralisedParameters)
                if getfield(a.p, kp) != getfield(b.p, kp)
                    return false
                end    
            end
        else
            if getfield(a, k) != getfield(b, k)
                return false
            end 
        end
    end
    return true
end

"""
    check_web!(com)

Remove double and canabalistic links from a community web
"""
function check_web!(com::Community)
    A = com.A
    sp_vec = com.sp
    #remove double links
    double_links = findall(UpperTriangular(A) .* LowerTriangular(A)' .!= 0)
    #loop over and clean double links, larger consumer stays
    for link = double_links
        i,j = link.I
        if sp_vec[i].n > sp_vec[j].n
            A[j,i] = 0
        else
            A[i,j] = 0
        end
    end

    #remove canabalism
    A[diagind(A)] .= 0

    #to do
    #remove isolated nodes    
    #ensure producer still exists!
end


"""
    community(sp_vec::Vector{Species}; T::Float64 = 0.5, R::Float64 = 42.0)

Generates an adjacency matrix for a given set of species using the niche model. removes spare nodes...
"""
function community(sp_vec::Vector{Species}; T::Float64 = 0.5, R::Float64 = 42.0)
    #add single producer
    # if !any(1e-2 .== [s.n for s = sp_vec])
    #     push!(sp_vec, species(0.0,T,n = 1e-2))
    # end


    N = length(sp_vec)
    A = zeros(N,N)
    n = zeros(N)

    #loop over species
    for (i,sp_i) = enumerate(sp_vec)
        n[i] = sp_i.n
        for (j,sp_j) = enumerate(sp_vec)
            #if j is within range of i
            if sp_j.n < (sp_i.c + (sp_i.r/2))
                if  sp_j.n > (sp_i.c - (sp_i.r/2))
                        A[i,j] = 1
                end
            end
        end
    end

    #ids
    ids = [x.id for x = sp_vec]

    com = Community(N, A, sp_vec, ids, n, T, R)
    check_web!(com)
    return com
end

"""
    community(sp_vec::Vector{Species}, N::Int64; T::Float64=0.5, T_range::Float64=0.1, R::Float64=42.0)

Randomly assembles a food web from a set of species of size N at temperature T. Species are selected within a range of T_range. Returns adjacency matrix. 
"""
# function community(sp_vec::Vector{Species}, N::Int64; T::Float64=0.5, T_range::Float64=0.1, R::Float64=42.0)
#     #select sp based on Temp range...
#     indx = findall([s.Tpk < (T + T_range) && s.Tpk > T - T_range for s = sp_vec])
#     indx = rand(indx, min(length(indx),N))
#     # N_T = 


#     #sample
#     # sp_vec_indx = rand(sp_vec_temp, N_T)

#     return community(sp_vec[indx], T=T, R=R)
# end

"""
    community(N,C)

Constructs a food web of size N and connectance C by randomly generating species. This function will automatically remove isolated species and add new ones. 
"""


"""
    parameterised_community(com::Community,p::GeneralisedParameters)

Adds generalised parameter set to a community
"""
function parameterised_community(com::Community,p::GeneralisedParameters)
    ParameterisedCommunity(com.N, com.A,com.sp,com.ids,com.n, com.T,com.R, p)
end

"""
    community(com::ParameterisedCommunity)

Coverts ParameterisedCommunity in Community object.
"""
function community(com::ParameterisedCommunity)
    Community(com.N, com.A,com.sp,com.ids,com.T,com.R)
end

#Community generation
"""
    add_species(com::Community, sp::Species)

Add species `sp` to a community and partially construct the niche web. 
"""
function add_species(com::Community, sp::Species)
    # get new sp list
    sp_vec_new = vcat(com.sp, [sp])
    n_new = vcat(com.n, [sp.n])

    #ids
    ids = [x.id for x = sp_vec_new]

    #add sp without recalculating web structure
    #new array
    A = zeros(size(com.A) .+ 1) 
    A[1:end-1 , 1:end-1] .= com.A
    #loop over sp
    for i = 1:size(com.A)[1] 
        #consumption
        a = 0
        if (com.sp[i].n < sp.c + (sp.r/2) && com.sp[i].n > sp.c - (sp.r/2))
            A[end,i] = 1
        end

        #predation
        if (sp.n < com.sp[i].c + (com.sp[i].r/2) && sp.n > com.sp[i].c - (com.sp[i].r/2))
            A[i,end] = 1
        end
    end

    new_com = Community(com.N + 1 , A, sp_vec_new, ids, n_new, com.T, com.R)
    # check_web!(new_com)

    return new_com
end

"""
    add_species(com::ParameterisedCommunity, sp::ParameterisedSpecies)

Adds a ParameterisedSpecies to a ParameterisedCommunity
"""
function add_species(com::ParameterisedCommunity, p_sp::ParameterisedSpecies)
    #add species
    com_new = add_species(community(com), Species(p_sp.n,p_sp.r,p_sp.c,p_sp.Tpk,p_sp.id))

    #calcualte new params
    N = size(com_new.A)[1]
    n = vcat(com.p.n, p_sp.n)
    α = vcat(com.p.α, (com.R .^ p_sp.n) .^ (-0.25))

    #update struct params
    sp = structural_parameters(com_new.A)
    
    #add exponential parameters
    ep = get_exponential_params(com.p)
    γ = vcat(ep.γ, p_sp.p.γ)
    λ = ones(N,N)
    μ = vcat(ep.μ, p_sp.p.μ)
    ϕ = vcat(ep.ϕ, p_sp.p.ϕ)
    ψ = vcat(ep.ψ, p_sp.p.ψ)

    #generate new paramstruct
    p_new = GeneralisedParameters(sp.N,n,com_new.A,α,sp.β,sp.χ,sp.ρ,sp.ρ̃,sp.σ,sp.σ̃,γ,λ,μ,ϕ,ψ)

    return parameterised_community(com_new, p_new)
end

"""
    remove_species(com::Community, id::UUID)

Remove species identified with `id`.
"""
function remove_species(com::Community, id::UUID)
    @assert id in com.ids "id not in community"

    indx = com.ids .!= id
    n_sp = com.sp[indx]
    n_ids = com.ids[indx]
    n_n = com.n[indx]

    new_com = Community(com.N - 1, com.A[indx,indx],n_sp,n_ids,n_n,com.T, com.R)
    # check_web!(new_com)

    return new_com
end

"""
    remove_species(com::ParameterisedCommunity, id::UUID)

Remove a species from a ParameterisedCommunity and update in the parameter set. 
"""
function remove_species(com::ParameterisedCommunity, id::UUID)
    @assert id in com.ids "id not in community"
    indx = com.ids .!= id

    com_new = remove_species(community(com),id)

     #calculate new params
     N = size(com_new.A)[1]
     n = com.p.n[indx]
     α = com.p.α[indx]

    sp = structural_parameters(com_new.A)

    #remove exponential parameters
    ep = get_exponential_params(com.p)

    γ = ep.γ[indx]
    λ = ones(N,N)
    μ = ep.μ[indx]
    ϕ = ep.ϕ[indx]
    ψ = ep.ψ[indx]
    
    #generate new paramstruct
    p_new = GeneralisedParameters(N,n,com_new.A,α,sp.β,sp.χ,sp.ρ,sp.ρ̃,sp.σ,sp.σ̃,γ,λ,μ,ϕ,ψ)
    
    return parameterised_community(com_new, p_new)
end

# """
#     remove_species(com::Community, sp::Species)

# Remove species `sp`
# """
# function remove_species(com::Community, sp::Species)
#     return remove_species(com, sp.id)
# end

"""
    move_species(com1::Community,com2::Community,id::UUID)

Move species identified by `id` from `com1` to `com2`
"""
function move_species(com1::Community,com2::Community, id::UUID)
    @assert id in com1.ids "no species with id in com1"
    @assert !(id in com2.ids) "id aleady in com2"
    #add to community 2
    com2_new = add_species(com2, com1.sp[findfirst(com1.ids .== id)])

    #remove from community 1
    com1_new = remove_species(com1, id)

    return com1_new, com2_new
end

function move_species(com1::ParameterisedCommunity,com2::ParameterisedCommunity, id::UUID)
    @assert id in com1.ids "no species with id in com1"

    #create parameterised_species
    p_sp = parameterised_species(com1, id)

    #add to community 2
    com2_new = add_species(com2, p_sp)

    #remove from community 1
    com1_new = remove_species(com1, id)

    return com1_new, com2_new
end

