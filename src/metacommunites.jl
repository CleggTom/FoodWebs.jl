"""
    MetaCommunity

type containing parameters defining a metacommunity including the communties `coms` as well as the matrix `D` which defines "distance" between communties.
"""
struct MetaCommunity
    coms::Array{Community}
    D::Array{Float64}
    T_mat::Array{Float64}
    sp::Vector{Species}
    sp_id::Vector{UUID}
    sp_loc::Dict{UUID,Vector{Int}}
    R::Float64
end

Base.show(io::IO, mc::MetaCommunity) =  print(io, "MetaCommunity M:", length(mc.coms)," Sp: ", length(mc.sp))

#constructers
"""
    metacommuntiy(coms::Array{Community})

Creates a MetaCommunity object from communities in `coms`.
"""
function metacommuntiy(coms::Array{Community})
    #get sp_vector
    sp_vec = unique(vcat([x.sp for x = coms]...))
    T_mat = [c.T for c = coms]

    #generate Sp location dictionaries
    sp_dict = Dict{UUID, Vector{Int}}()
    for i = eachindex(coms)
        for id = coms[i].ids
            if !haskey(sp_dict, id)
                sp_dict[id] = [i]
            else
                push!(sp_dict[id], i)
            end
        end
    end

    #generate temperature distance matrix
    D = zeros(length(T_mat), length(T_mat))
    for i = eachindex(T_mat)
        for j = eachindex(T_mat)
            #calculate distance matrix
            D[i,j] = T_mat[j] - T_mat[i]
            D[i,j] = D[i,j] <= 0 ? Inf : D[i,j]
        end
    end

    #get sp vecs
    id_vec_mc = collect(keys(sp_dict))
    sp_vec_mc = hcat([sp_vec[findall([id == s.id for s = sp_vec])] for id = id_vec_mc]...)[:]

    #get communtiy Ratios
    x = [com.R for com = coms]
    @assert all( y -> y==x[1], x) "ppmr must be the same across communities"

    return MetaCommunity(coms, D, T_mat, sp_vec_mc, id_vec_mc, sp_dict, coms[1].R)
end



"""
    metacommuntiy(sp_vec::Vector{Species}, N::Int64, T_mat, T_range::Float64, R::Float64)

Generate metacommuntiy from a set of species with community size N and temperatures T_mat. Sampling is done with species T_range unsits of T. 
"""
function metacommuntiy(sp_vec::Vector{Species}, N::Int64, T_mat; T_range::Float64=0.1, R::Float64=42.0)
    coms = [community(sp_vec, N, T=T, T_range=T_range, R=R) for T = T_mat]

    return metacommuntiy(coms)
end


"""
    stable_metacommunity(sp_vec::Vector{Species}, N::Int64, T_mat, T_range::Float64, R::Float64, psw_threshold::Float64 = 0.9, N_trials::Int = 100, max_draws::Int = 100)

Contstruct a stable metacommunity by resampling communties till they are stable. The metacommunity is generated from sp_vec. 
"""
function stable_metacommunity(sp_vec::Vector{Species}, N::Int64, T_mat; T_range::Float64 = 0.1, R::Float64=42.0, psw_threshold::Float64 = 0.9, N_trials::Int = 100, max_draws::Int = 100, verbose = false, vk = 5)
    #inital communty generation
    mc = metacommuntiy(sp_vec, N, T_mat, T_range=T_range, R=R)
    psw = proportion_stable_webs(mc , N_trials = N_trials)
    coms = mc.coms

    #resample unstable communties
    k = 0
    while (k < max_draws) && any(psw .< psw_threshold)
       #replace unstable communtiy
        Threads.@threads for i = eachindex(psw)
            if psw[i] .< psw_threshold
                coms[i] = community(sp_vec, N , T=coms[i].T, T_range=T_range, R=R) 
                psw[i] = proportion_stable_webs(coms[i], N_trials = N_trials)
            end
        end
        
        verbose && k % vk == 0 && println("draw:", k, " psw: ", psw)

        k += 1
    end

    mc = metacommuntiy(coms)

    return(mc)
end

#functions
"""
    move_sp_meta!(mc::MetaCommunity, a, b, id)

Moves species with `id` from community a to b in a metacommuntiy. a and b are the  1-D indexes of the communties. 
"""
function move_sp_meta!(mc::MetaCommunity, a, b, id)
    @assert !(id in mc.coms[b].ids) "Species with `id` is already in communtiy b"
    mc.coms[a], mc.coms[b] = move_species(mc.coms[a], mc.coms[b], id)

    #update loc
    append!(mc.sp_loc[id], b)
    filter!(x -> x != a, mc.sp_loc[id])
end

#single random movement
"""
    random_dispersal!(mc, p_dispersal = :weighted, d_dispersal = :weighted)

Randomly selects and moves a single Species between webs in a MetaCommunity. The movement is split into two stages:

1) Select Species
Select the species to move based on relative body size `n`. This is done by sampling with weighted probablities set to 1 - exp(-n ^ 0.75)

2) Select site
Select the site the species will disperse to. This is done by considering the distance matrix 
"""
# function random_dispersal!(mc; p_dispersal = :p, d_dispersal = :p)
#     @assert p_dispersal ∈ [:p, :r]
#     @assert d_dispersal ∈ [:p, :r]

#     #sample sp to disperse
#     if p_dispersal == :p
#         id = sample(mc.sp_id,  Weights([1 - exp(-sp.n ^ 0.75) for sp = mc.sp]))
#     elseif p_dispersal == :r
#         id = sample(mc.sp_id)
#     end

#     #get from location
#     from = sample(mc.sp_loc[id])

#     if d_dispersal == :p
#         #get to
#         #distance rate
#         n = mc.sp[findall(id .== mc.sp_id)[1]].n
#         λ = (1 - n)^(0.75)
#         w = exp.(-λ * mc.D[from,:] * 2)

#     elseif d_dispersal == :r
#         w = .!isinf.(mc.D[from,:])
#     end

#     #skip if nowhere to disperse to
#     if all(w .== 0.0)
#         # print("cant disperse")
#         return 0,0
#     end

#     to = sample(1:length(mc.T_mat), Weights(w))

#     if !(id in mc.coms[to].ids)
#         move_sp_meta!(mc, from, to, id)
#     end
    
#     return from, to
# end


"""
    multiple_dispersal!(mc; p_dispersal = :p, d_dispersal = :p, K = 5)

Disperse multiple species at once. Modifies a MetaCommunity object in place. The process is split into two stages:
1) Species selection

2) Communtiy selection
"""
function multiple_dispersal!(mc; p_dispersal = :p, d_dispersal = :p, K = 5, verbose = false)
    @assert p_dispersal ∈ [:k, :p, :r]
    @assert d_dispersal ∈ [:p, :r]
    check_metacommunity(mc)

    #calculate mass + rate
   
    M = [mc.R ^ sp.n for sp = mc.sp]
    λ = M .^ 0.75

    #sample species to disperse
    if p_dispersal == :r #uniform random
        ids = sample(mc.sp_id, K, replace = false)
    elseif p_dispersal == :k #K samples
        ids = sample(mc.sp_id,  Weights(1 .- exp.( -λ)), K)
    elseif p_dispersal == :p #probabalistic
        p = 1 .- exp.( -λ)
        ids = mc.sp_id[findall(rand(length(p)) .< p)]
    end

    #get from locations
    from = sample.(get.(Ref(mc.sp_loc), ids, 0))
    
    if verbose 
        println(" moving:", length(from))
    end
    
    #sample distances
    for k = eachindex(from) #loop over species
        if d_dispersal == :r
            #randomly sample
            to = sample(1: length(mc.D[from[k],:]))
        elseif d_dispersal == :p
            #sample from exponential Distribution upwards
            d = rand(truncated(Exponential(λ[k]), 0, 1))
            to = findmin((mc.D[from[k],:] .- d).^2)[2]
        end
        #if sp is not there add it
        if !(ids[k] in mc.coms[to].ids)
            move_sp_meta!(mc, from[k], to, ids[k])
        end
        check_metacommunity(mc)
    end

end

# # #1) GENERATE SPECIES POOL
# N_pool = 100000
# sp_vec = [species(0.1) for i = 1:N_pool];

# #2) ASSEMBLE COMMUNITIES
# #set temperatures
# N_T = 2
# t_vec = range(0,1,length = N_T)

# #generate metacommunty
# mc1 = stable_metacommunity(sp_vec, 10, t_vec, T_range = 0.1, R = 43.0,
#          psw_threshold = 0.9, max_draws = 100, verbose = true, vk = 1)

# mc2 = deepcopy(mc1)
# check_metacommunity(mc1)

# id = mc1.sp_id[1]
# from = mc1.sp_loc[mc1.sp_id[1]][1]
# to = [1,2][findfirst(from .∉ [1, 2])]

# move_sp_meta!(mc1, from, to, id)

# check_metacommunity(mc1)
# check_metacommunity(mc2)

# multiple_dispersal!(mc1)

#move


"""
    test_metacommunity(mc)

Tests to see if the sp_loc dictionary is correct (i.e. all species are where they should be).
"""
function check_metacommunity(mc::MetaCommunity)
    for (i,sp) in enumerate(mc.sp)
        coms = mc.sp_loc[sp.id]
        for c in coms
            @assert sp in mc.coms[c].sp "$sp is not in $c"
        end
    end
end

"""
    get_richness(mc::MetaCommunity)

Gets vector of richness in a metacommuntiy
"""
function get_richness(mc::MetaCommunity)
    [length(c.sp) for c = mc.coms]
end