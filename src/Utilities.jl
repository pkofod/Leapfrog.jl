K(c::Float64, cost::CostParameters) = cost.k1/(1+cost.k2*c)
Φ(vN, vI) = maximum([vN vI])
r(ci::Float64, cj::Float64) = ci >= cj ? 0.0 : cj - ci
r1(c1, c2) = r(c1::Float64, c2::Float64)
r2(c1, c2) = r(c2::Float64, c1::Float64)
number_of_strings(nC::Int64) = div(nC*(nC+1)*(2*nC+1),6)

# Functions for the lowest layer

H1(ic1, ic2, pi, τ, ess)     = pi*Φ(τ[ess.nC].EQs[ic1, ic2, ess[ic1,ic2,ess.nC]].vN1,τ[ess.nC].EQs[ic1, ic2, ess[ic1,ic2,ess.nC]].vI1)
H2(ic1, ic2, pi, τ, ess)     = pi*Φ(τ[ess.nC].EQs[ic1, ic2, ess[ic1,ic2,ess.nC]].vN2,τ[ess.nC].EQs[ic1, ic2, ess[ic1,ic2,ess.nC]].vI2)

# Functions for the upper layers

H1(ic1, ic2, ic, pi, τ, ess) = (1-pi)*Φ(τ[ic+1].EQs[ic1, ic2, ess[ic1,ic2,ic+1]].vN1,τ[ic+1].EQs[ic1, ic2, ess[ic1,ic2,ic+1]].vI1) +
                               pi*Φ(τ[ic].EQs[ic1, ic2, ess[ic1,ic2,ic]].vN1,τ[ic].EQs[ic1, ic2, ess[ic1,ic2,ic]].vI1)
H2(ic1, ic2, ic, pi, τ, ess) = (1-pi)*Φ(τ[ic+1].EQs[ic1, ic2, ess[ic1,ic2,ic+1]].vN2,τ[ic+1].EQs[ic1, ic2, ess[ic1,ic2,ic+1]].vI2) +
                               pi*Φ(τ[ic].EQs[ic1, ic2, ess[ic1,ic2,ic]].vN2,τ[ic].EQs[ic1, ic2, ess[ic1,ic2,ic]].vI2)
# Quadratic solver

function quad(a::Float64, b::Float64, c::Float64)
    d = b^2 - 4*a*c
    if abs(a) < 1e-8
        pstar = [0. ; 1.;-c/b]
    else
        if d < 0
            pstar = [0. ;1.]
        elseif d == 0.
            pstar = [0. ;1. ; -b/2a]
        else
            pstar = [0. ;1. ; (-b - sqrt(d))/(2a); (-b + sqrt(d))/(2a)]
        end
    end
    return pstar
end

function jump!(x::ESS)
    activated = false
    if x[1,1,1] == -1
        throw(error("This ESS has already overflown!"))
    end
    X = Array(Int64,length(x))
    R = 1 # Remainder; first "remainder" is one, since we add one
    @inbounds for i in reverse(1:length(x)) # Slightly inefficient, as we go through R = 0 cases
        X[i] = mod(x.number[i] + R,x.bases[i])
        R = div(x.number[i]+R,x.bases[i])
        if R == 0 && activated == false
            x.active = i
            activated = true
        end
    end
    if R > 0
        println("No more equilibria to check.")
        x.number = -1*ones(Int64,length(x))
    else
        x.number = X
    end
end
