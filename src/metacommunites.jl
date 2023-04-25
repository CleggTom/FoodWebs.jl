
#constructers
"""
    metacommuntiy(coms::Array{Community})

Creates a MetaCommunity object from communities in `coms`.
"""
function metacommuntiy(coms::Array{T}) where T <: AbstractCommunity
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
        for i = eachindex(psw)
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

function stable_metacommunity_p(sp_vec::Vector{Species}, N::Int64, T_mat; T_range::Float64 = 0.1, R::Float64=42.0, N_trials::Int = 100, max_draws::Int = 100, verbose = false, vk = 5)
        #inital communty generation
        mc = metacommuntiy(sp_vec, N, T_mat, T_range=T_range, R=R)
        coms = Vector{AbstractCommunity}(undef, length(T_mat))

        for (i,c) = enumerate(mc.coms)
            coms[i] = stable_parameterisation(c, N_trials)
        end
        
        #resample unstable
        k = 0
        while (k < max_draws) && any(isa.(coms,Ref(Community)))

            for (i,c) = enumerate(coms)
                if isa(c, ParameterisedCommunity)
                    continue
                end
                com_new = community(sp_vec, N , T=coms[i].T, T_range=T_range, R=R) 
                coms[i] = stable_parameterisation(com_new, N_trials)
            
            end

            verbose && k % vk == 0 && println("draw:", k, " stable: ", typeof.(coms) .== ParameterisedCommunity)

            k += 1
        end

        mc_p = metacommuntiy(coms)

        return(mc_p)
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

    #calculate mass +ß rate
   
    M = [mc.R ^ sp.n for sp = mc.sp]
    d0 = -log(1 - 0.9) / (mc.R^0.75)
    λ_p = d0 .* (M .^ 0.75)

    #sample species to disperse
    if p_dispersal == :r #uniform random
        ids = sample(mc.sp_id, K, replace = false)
    elseif p_dispersal == :k #K samples
        ids = sample(mc.sp_id,  Weights(1 .- exp.( -λ_p)), K)
    elseif p_dispersal == :p #probabalistic
        p = 1 .- exp.( -λ_p)
        ids = mc.sp_id[findall(rand(length(p)) .< p)]
    end

    #get from locations
    from = sample.(get.(Ref(mc.sp_loc), ids, 0))
    
    if verbose 
        println(" moving:", length(from))
    end
    
    #sample distances
    d0 = -log(1 - 0.9) / 0.5(mc.R^-0.75)
    λ_d = d0 .* (M .^ -0.75)
    for k = eachindex(from) #loop over species
        if d_dispersal == :r
            #randomly sample
            to = sample(1: length(mc.D[from[k],:]))
        elseif d_dispersal == :p
            #sample from Exponential Distribution upwards
            d = rand(Exponential(1 / λ_d[k]))
            to = findmin((mc.D[from[k],:] .- d).^2)[2]
        end
        #if sp is not there add it
        if !(ids[k] in mc.coms[to].ids)
            move_sp_meta!(mc, from[k], to, ids[k])
        end
        check_metacommunity(mc)
    end

end

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