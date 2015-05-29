module Leapfrog

import Base.length
export init_cost, init, init_all, solve, getindex


include("Types.jl")
include("Utilities.jl")
include("Last.jl")
include("Upper.jl")

function solve(cost, mp)
    ess = ESS(cost.nC) # Start with feasible ESS; we're solving for all equilibria
    number_of_equilibria=0
    τ = Array(Stage, cost.nC)
    for i in reverse(1:cost.nC)
        τ[i] = Stage(Array(EQ, i, i, 5))
    end

    solve_last_corner!(cost, mp, τ)
    solve_last_edge!(cost, mp, τ)
    solve_last_interior!(cost, mp, τ, ess)
    for ic = reverse(1:cost.nC-1)
        solve_corner!(cost, mp, ic, τ, ess)
        if ic != cost.nC
            solve_edges!(cost, mp, ic, τ, ess)
            solve_interior!(cost, mp, ic, τ, ess)
        end
    end

    number_of_equilibria += 1
    jump!(ess)
    if ess[1,1,1] != 0
        for void=1:100000
            
            if ess.ess_to_ic[ess.active][3] == ess.nC
                solve_last_interior!(cost, mp, ess.ess_to_ic[ess.active][1], ess.ess_to_ic[ess.active][2], τ, ess)
            elseif ess.ess_to_ic[ess.active][3] < ess.nC
                solve_interior!(cost, mp, ess.ess_to_ic[ess.active][3], ess.ess_to_ic[ess.active][1], ess.ess_to_ic[ess.active][2], τ, ess)
            end 

            for ic = reverse(1:ess.ess_to_ic[ess.active][3]-1) # already solved stage!
                solve_corner!(cost, mp, ic, τ, ess)
                if ic != cost.nC
                    solve_edges!(cost, mp, ic, τ, ess)             
                    solve_interior!(cost, mp, ic, τ, ess)
                end
            end
            number_of_equilibria+=1
            jump!(ess)
            if ess[1,1,1] == 0
                break
            end
        end
    end
    return τ, ess, number_of_equilibria
end

function init_cost(switch::Int64)
    if switch == 1
        k1 = 8.3
        k2 = 1.
        C  = linspace(5, 0, 4)
        Cmin = minimum(C)
        nC = length(C)
        p = [1.; 1.; 1.; 0]
        π = 1-p
        return CostParameters(k1, k2, C, Cmin, nC, p, π)
    else
        k1 = 60.3
        k2 = 1.
        C  = linspace(5, 0, 3)
        Cmin = minimum(C)
        nC = length(C)
        p = [1.; 1.; 0.]
        π = 1-p
        return CostParameters(k1, k2, C, Cmin, nC, p, π)
    end
end

function init_mp(cost)
    R = 0.05
    Δt = 1
    β = exp(- R*Δt) # R is continuous time iterest rate
    T = (cost.nC-1)*3 + 1
    return RLSModel(R, β, Δt, T, number_of_strings(cost.nC))
end

function init_all(switch::Int64)
    cost = init_cost(switch)
    return init_mp(cost), cost
end
end