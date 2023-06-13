
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
            D[i,j] = D[i,j] < 0 ? Inf : D[i,j]
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
    #create temperature lookup
    T_vals = unique(T_mat) 
    T_lookup = Dict(T => Species[] for T = T_vals)
    for sp = sp_vec
        for T = T_vals
            if sp.Tpk < (T + T_range) && sp.Tpk > T - T_range
                push!(T_lookup[T], sp)
            end
        end
    end


    #proportion stable webs
    psw = zeros(size(T_mat))
    coms = similar(T_mat, Community)
    J = zeros(N,N)

    #resample communties
    k = 0
    while (k < max_draws) && any(psw .< psw_threshold)
       #replace unstable communtiy
        for i = eachindex(T_mat)
            if psw[i] .< psw_threshold
                sp_indx = rand(T_lookup[T_mat[i]], N)
                coms[i] = community(sp_indx, T = T_mat[i], R = R) 
                psw[i] = proportion_stable_webs!(J, coms[i], N_trials = N_trials)
            end
        end
        
        verbose && k % vk == 0 && println("draw:", k, " psw: ", psw)

        k += 1
    end

    mc = metacommuntiy(coms)

    return(mc)
end

function stable_metacommunity(N::Int64, C::Float64, T_mat; T_range::Float64 = 0.1, R::Float64=42.0, psw_threshold::Float64 = 0.9, N_trials::Int = 100, max_draws::Int = 100, verbose = false, vk = 5)
    #proportion stable webs
    psw = zeros(size(T_mat))
    coms = similar(T_mat, Community)
    J = zeros(N,N)

     #resample communties
     k = 0
     while (k < max_draws) && any(psw .< psw_threshold)
        #replace unstable communtiy
         for i = eachindex(T_mat)
             if psw[i] .< psw_threshold
                 sp_indx = species(C, fill(T_mat[i], N), N)
                 coms[i] = community(sp_indx, T = T_mat[i], R = R) 
                 #remove isolated 
                #  isolated = 
                 psw[i] = proportion_stable_webs!(J, coms[i], N_trials = N_trials)
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

function remove_sp_meta!(mc::MetaCommunity, a, id)
    mc.coms[a] = remove_species(mc.coms[a], id)
    
    #update loc
    filter!(x -> x != a, mc.sp_loc[id])
end

"""
    move_sp_meta!(mc::MetaCommunity, a, b, id)

Moves species with `id` from community a to b in a metacommuntiy. a and b are the  1-D indexes of the communties. 
"""
function move_sp_meta!(mc::MetaCommunity, a, b, id)
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
function multiple_dispersal!(mc; p_dispersal = :p, d_dispersal = :p, verbose = false, αp = 0.75, αd = 0.75)
    @assert p_dispersal ∈ [:p, :r, :a]
    @assert d_dispersal ∈ [:p, :r]
    check_metacommunity(mc)


    to_disperse = Vector{Tuple{UUID, Pair{Int,Int}}}()
    #Loop over communities
    for (f,c) = enumerate(mc.coms)
        M = c.R .^ c.n #calculate mass
        
        #calculate probability of dispersal in C
        if p_dispersal == :a
            p = ones(c.N)
        elseif p_dispersal == :r
            p = fill(mean(c.n), c.N)
        elseif p_dispersal == :p
            # λ = 0.05 .* M .^ αp
            # p = 1 - exp(-λ)
            p = c.n
        end

        #test who will disperse
        to_d = rand(c.N) .< p

        for sp = 1:c.N
            if to_d[sp]
                #calculate dispersal distance
                if d_dispersal == :p
                    d = 0.5 * (M[sp]/mc.R)^αd
                    t = findmin(abs.(mc.D[f,:] .- d))[2]
                elseif d_dispersal == :r
                    t = rand(f:length(mc.T_mat))
                end

                push!(to_disperse, (c.ids[sp], f => t))
            end
        end
    
    end


    for td = to_disperse
        f,t,id = td[2][1], td[2][2], td[1]
        #if sp is not there add it
        if t != length(mc.T_mat)
            if !(id in mc.coms[t].ids) && (id in mc.coms[f].ids)
                move_sp_meta!(mc, f, t, id)
            end
        else
            remove_sp_meta!(mc,f,id)
        end
    end
    
    [check_web!(c) for c = mc.coms]
    check_metacommunity(mc)

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