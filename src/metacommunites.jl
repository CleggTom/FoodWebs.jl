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
end

Base.show(io::IO, mc::MetaCommunity) =  print(io, "MetaCommunity M:", length(mc.coms)," Sp: ", length(mc.sp))

#constructers
"""
    metacommuntiy(sp_vec::Vector{Species}, N::Int64, T_mat, T_range::Float64)

Creates a MetaCommunity object from species in `sp_vec` with each community of size N. Community temperatures are defined by the `T_mat` argument. 
"""
function metacommuntiy(sp_vec::Vector{Species}, N::Int64, T_mat, T_range::Float64)
    #generate communties
    coms = [community(sp_vec, N, T, T_range) for T = T_mat]
    
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

    D = zeros(length(T_mat), length(T_mat))
    for i = 1:length(T_mat)
        for j = 1:length(T_mat)
            #calculate distance matrix
            D[i,j] = T_mat[j] - T_mat[i]
            D[i,j] = D[i,j] <= 0 ? Inf : D[i,j]
        end
    end

    #get sp vecs
    id_vec_mc = collect(keys(sp_dict))
    sp_vec_mc = hcat([sp_vec[findall([id == s.id for s = sp_vec])] for id = id_vec_mc]...)[:]


    return MetaCommunity(coms, D, T_mat, sp_vec_mc, id_vec_mc, sp_dict)
end

function stable_metacommunity(sp_vec::Vector{Species}, N::Int64, T_mat, T_range::Float64, psw_threshold::Float64 = 0.9)
    mc = metacommuntiy(sp_vec, N, T_mat, T_range)
    psw = proportion_stable_webs(mc ,N = 100)

    while any(psw .< psw_threshold)
        mc = metacommuntiy(sp_vec, N, T_mat, T_range)
        psw = proportion_stable_webs(mc,N = 100)
    end

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
    random_dispersal(mc)

Randomly selects and moves a single Species between webs in a MetaCommunity. The movement is split into two stages:

1) Select Species
Select the species to move based on relative body size `n`. This is done by sampling with weighted probablities set to 

"""
function random_dispersal(mc)
    #sample sp to disperse
    id = sample(mc.sp_id,  Weights([1 - exp(-sp.n ^ 0.75) for sp = mc.sp]))

    #get from location
    from = sample(mc.sp_loc[id])

    #get to
    #distance rate
    n = mc.sp[findall(id .== mc.sp_id)[1]].n
    λ = (1 - n)^(0.75)
    w = exp.(-λ * mc.D[from,:] * 2)

    if all(w .== 0.0)
        print("cant disperse")
        return
    end

    to = sample(1:length(mc.T_mat), Weights(w))
    print("moving sp: ",id, " from ",from, " to ", to)
    @assert from < to "Oh No"
    if !(id in mc.coms[to].ids)
        move_sp_meta!(mc, from, to, id)
    end
    
end

sp_vec = [species(0.1) for i = 1:5000]
T_mat = range(0,1, length = 10)
T_range = 0.1
N = 100

mc = metacommuntiy(sp_vec, N, T_mat, T_range)

@time proportion_stable_webs(mc, 10)
#check sp locations

function check_metacommunity(mc)
    for (i,sp) in enumerate(mc.sp)
        coms = mc.sp_loc[sp.id]
        for c in coms
            @assert sp in mc.coms[c].sp "$c $sp"
        end
    end
end