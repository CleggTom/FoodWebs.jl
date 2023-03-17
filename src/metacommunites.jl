
"""
    MetaCommunity

type containing parameters defining a metacommunity including its dispersal matrix B.
"""
struct MetaCommunity
    coms::Vector{Community}
    B::Matrix{Float64}
end

