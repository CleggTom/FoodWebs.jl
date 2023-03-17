"""
    Species

Type to store niche model parameters for a single species, its thermal optima and a unique ID.  
"""
struct Species
    n::Float64
    r::Float64
    c::Float64
    Tpk::Float64
    id::UUID
end

"""
    species(C::Float64)

Constructor function for Sp type. Randomly generates niche parameters given a connectance value
"""
function species(C::Float64)
    #niche params
    β = (1 / (2*C)) - 1
    n = rand(Uniform(0,1))
    r = n * rand(Beta(1.0, β))
    c = rand(Uniform(r / 2, n))
    #extra params
    Tpk = rand()
    uuid = UUIDs.uuid1()
    return(Species(n,r,c,Tpk,uuid))
end

"""
    Community

Type containing parameters defining a community including its interaction matrix and a vector of species within.
"""
struct Community
    A::Matrix{Float64}
    sp::Vector{Species}
    ids::Vector{UUID}
    T::Float64
end

"""
    isequalcom(a::Community, b::Community)

Test if two community objects are equal. This forces comparison by equality (`==`) as opposed to identity ('===') which the default `isequal` method uses. 
"""
function isequalcommunity(a::Community, b::Community)
    #get keys
    com_keys = fieldnames(Community)
    #test if all fields are the same
    return getfield.(Ref(a), com_keys) == getfield.(Ref(b), com_keys) 
end

Base.show(io::IO, com::Community) = print(io, "Community with ", length(com.sp)," species")

"""
    community(sp_vec::Vector{Sp}, T::Float64)

Generates an adjacency matrix for a given set of species using the niche model. 
"""
function community(sp_vec::Vector{Species}, T::Float64)
    N = length(sp_vec)
    A = zeros(N,N)

    for (i,sp_i) = enumerate(sp_vec)
        for (j,sp_j) = enumerate(sp_vec)
            if sp_j.n < (sp_i.c + (sp_i.r/2))
                if  sp_j.n > (sp_i.c - (sp_i.r/2))
                    if i ≠ j
                        A[i,j] = 1
                    end
                end
            end
        end
    end

    ids = [x.id for x = sp_vec]

    return(Community(A, sp_vec, ids, T))
end

"""
    community(sp_vec::Vector{Sp}, N::Int64, T::Float64, T_range::Float64)

Randomly assembles a food web from a set of species of size N at temperature T. Species are selected within a range of T_range. Returns adjacency matrix. 
"""
function community(sp_vec::Vector{Species}, N::Int64, T::Float64, T_range::Float64)
    #select sp based on Temp range...
    indx = findall([s.Tpk < (T + T_range) && s.Tpk > T - T_range for s = sp_vec])
    
    N_T = min(length(indx),N)

    sp_vec_temp = sp_vec[indx]

    #sample
    sp_vec_indx = sp_vec_temp[sortperm(rand(length(indx)))[1:N_T]]

    community = community(sp_vec_indx)

    return(community)
end

"""
    community(N,C)

Constructs a food web of size N and connectance C by randomly generating species. This function will automatically remove isolated species and add new ones. 
"""

#dispersal


"""
    add_species(com::Community, sp::Species)

Add species `sp` to a community and recalculate the niche web. 
"""
function add_species(com::Community, sp::Species)
    # get new sp list
    sp_vec_new = vcat(com.sp, [sp])

    return(community(sp_vec_new, com.T))
end

"""
    remove_species(com::Community, id::UUID)

Remove species identified with `id`.
"""
function remove_species(com::Community, id::UUID)
    @assert id in com.ids

    indx = com.ids .!= id

    return(community(com.sp[indx], com.T))
end

"""
    remove_species(com::Community, sp::Species)

Remove species `sp`
"""
function remove_species(com::Community, sp::Species)
    return(remove_species(com, sp.id))
end

"""
    move_species(com1::Community,com2::Community,id::UUID)

Move species identified by `id` from `com1` to `com2`
"""
function move_species(com1::Community,com2::Community,id::UUID)
    #add to community 2
    com2_new = add_species(com2, com1.sp[findfirst(com1.ids .== id)])

    #remove from community 1
    com1_new = remove_species(com1, id)

    return(com1_new, com2_new)
end


"""
    move_species(com1::Community,com2::Community,sp::Species)

Move `sp`` from `com1` to `com2`. note that the species `sp` must be present in `com1`
"""
function move_species(com1::Community,com2::Community,sp::Species)
    @assert sp in com1.sp
    return(move_species(com1,com2,sp.id))
end