immutable RLSModel
    R::Float64
    β::Float64
    Δt::Float64
    T::Int64
    nESS::Int64
end

immutable CostParameters
    k1::Float64
    k2::Float64
    C::LinSpace # Capital C refers to the potential C, not a realization
    Cmin::Float64
    nC::Int64
    p::Vector{Float64}
    π::Vector{Float64}
end

type EQ
    P1::Float64
    vN1::Float64
    vI1::Float64
    P2::Float64
    vN2::Float64
    vI2::Float64
end

type Stage
    EQs::Array{EQ,3}
end

type ESS
    number::Array{Int64,1}
    bases::Array{Int64,1}
    nC::Int64
    nESS::Int64
    active::Int64
    ic_to_ess::Array{Int64,3}
    ess_to_ic::Dict{Int64, Array{Int64,1}}
end

function ESS(nC::Int64)
    a,b = index_dict(nC)
    ESS(zeros(Int64,number_of_strings(nC)), ones(Int64,number_of_strings(nC)),nC,number_of_strings(nC),1,a,b)
end


type Statistics
    P1::Float64
    P2::Float64
end

function getindex(x::ESS, ic1::Int64, ic2::Int64, ic::Int64)
    return x.number[x.ic_to_ess[ic1,ic2,ic]]+1
end

function index(nC::Int64, ic1::Int64, ic2::Int64, ic::Int64)
    if (ic1,ic2) == (ic,ic)
        return 1 + div(nC*(nC+1)*(2*nC+1),6) - div(ic*(ic+1)*(2*ic+1),6)
    elseif ic2 == ic
        return 1 + div(nC*(nC+1)*(2*nC+1),6) - div(ic*(ic+1)*(2*ic+1),6) + ic1
    elseif ic1 == ic
        return 1 + div(nC*(nC+1)*(2*nC+1),6) - div(ic*(ic+1)*(2*ic+1),6) + ic - 1 + ic2
    else
        return 1 + div(nC*(nC+1)*(2*nC+1),6) - div(ic*(ic+1)*(2*ic+1),6) + 2*(ic - 1) + sub2ind((ic-1,ic-1),ic1,ic2)
    end
end

function index_dict(nC::Int64)
    count = 1
    nESS = div(nC*(nC+1)*(2*nC+1),6)
    A = Dict{Int64, Array{Int64,1}}()
    B = Array{Int64}(nC, nC, nC)
    for ic = 1:nC
        for ic1 = 1:ic
            for ic2 = 1:ic
                A[index(nC, ic1, ic2, ic)] = [ic1,ic2,ic]
                B[ic1,ic2,ic] = index(nC, ic1, ic2, ic)
                count += 1
            end
        end
    end
    return B,A
end

length(x::ESS) = length(x.number)