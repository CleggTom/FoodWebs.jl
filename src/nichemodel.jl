"""
    Sp

type to store niche model parameters for a single species, its thermal optima and a unique ID.  
"""
struct Sp
    n::Float64
    r::Float64
    c::Float64
    Tpk::Float64
    id::UUID
end

"""
    Community

type containing parameters defining a community including its interaction matrix and a vector of species within.
"""
struct Community
    A::Matrix{Float64}
    sp::Vector{Sp}
end

"""
    sp(C::Float64)

Constructor function for Sp type. Randomly generates niche parameters given a connectance value
"""
function Sp(C::Float64)
    β = (1 / (2*C)) - 1
    n = rand(Uniform(0,1))
    r = n * rand(Beta(1.0, β))
    c = rand(Uniform(r / 2, n))
    Tpk = rand()
    uuid = UUIDs.uuid1()
    return(Sp(n,r,c,Tpk,uuid))
end


"""
    web(sp_vec::Vector{Sp})

Generates an adjacency matrix for a given set of species using the niche model. 
"""
function web(sp_vec::Vector{Sp})
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

    return(Community(A, sp_vec))
end

"""
    web(sp_vec::Vector{Sp}, N::Int64, T::Float64, T_range::Float64)

Randomly assembles a food web from a set of species of size N at temperature T. Species are selected within a range of T_range. Returns adjacency matrix. 
"""
function web(sp_vec::Vector{Sp}, N::Int64, T::Float64, T_range::Float64)
    #select sp based on Temp range...
    indx = findall([s.Tpk < T + T_range && s.Tpk > T - T_range for s = sp_vec])
    
    N_T = min(length(indx),N)

    #sample
    sp_vec_indx = sp_vec[sortperm(rand(length(indx)))[1:N_T]]

    community = web(sp_vec_indx)

    return(community)
end

"""
    web(N,C)

Constructs a food web of size N and connectance C by randomly generating species. This function will automatically remove isolated species and add new ones. 
"""




