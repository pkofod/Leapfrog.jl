function solve_corner!(cost::CostParameters, mp::RLSModel, ic::Int64, τ::Array{Stage,1}, ESS)
    c = cost.C[ic]
    P1 = 0.
    vN1 = mp.β*((1-cost.π[ic])*Φ(τ[ic+1].EQs[ic, ic, ESS[ic,ic,ic+1]].vN1, τ[ic+1].EQs[ic, ic, ESS[ic,ic,ic+1]].vI1)+Φ(0, -K(c, cost)))/(1 - cost.π[ic]*mp.β)
    vI1 = vN1 - K(c, cost)
    vN2 = mp.β*( (1-cost.π[ic])*Φ(τ[ic+1].EQs[ic, ic, ESS[ic,ic,ic+1]].vN2, τ[ic+1].EQs[ic, ic, ESS[ic,ic,ic+1]].vI2)+Φ(0, -K(c, cost)))/(1 - cost.π[ic]*mp.β)
    vI2 = vN2 - K(c, cost)
    P2 = 0.
    τ[ic].EQs[ic,ic,1] = EQ(P1, vN1, vI1, P2, vN2, vI2)
end

function solve_edges!(cost::CostParameters, mp::RLSModel, ic::Int64, τ::Array{Stage,1}, ess::ESS)
    c = cost.C[ic]
    pi = cost.π[ic]

    # if at (c1,c,c) edge
    for ic1 = 1:ic-1
        c1 = cost.C[ic1]
        vI1 = -K(c,cost) + mp.β*H1(ic,ic,ic, pi, τ, ess)

        #vN1search(x) = mp.β*((1-pi)*Φ(τ[ic+1].EQs[ic1,ic,ess[ic1,ic,ic+1]].vN1, τ[ic+1].EQs[ic1, ic, ess[ic1,ic,ic+1]].vI1)+pi*Φ(x,vI1)) - x
        #vN1 = fzero(vN1search,0;order=2)
        #vN1 = mp.β*vI1

        a = vI1
        A = mp.β*(1-pi)*Φ(τ[ic+1].EQs[ic1,ic,ess[ic1,ic,ic+1]].vN1, τ[ic+1].EQs[ic1, ic, ess[ic1,ic,ic+1]].vI1)
        B = mp.β*(1-pi)
        vN1 = 0
        if a == 0
            if  A >= 0
                vN1 = -A/(B-1)
            elseif A < 0
                vN1 = A
            end
        elseif a > 0
            if A >= a
                vN1 = -A/(B-1)
            elseif 0 < A < a
                if B <= (a-A)/a
                    vN1 = a*B+A
                elseif (a-A)/a < B
                    vN1 = -A/(B-1)
                end
            elseif A <= 0
                x = a*B+A
            end
        elseif a < 0
            if A >= 0
                vN1 = -A/(B-1)
            elseif a < A < 0
                if B <= (a-A)/a
                    vN1 = -A/(B-1)
                elseif (a-A)/a < B
                    vN1 = a*B+A
                end
            elseif A <= a
                vN1 = a*B+A
            end
        end

        P1 = 0.

        vN2 = (r2(c1,c) + mp.β*(P1*((1-pi)*Φ(τ[ic+1].EQs[ic, ic, ess[ic,ic,ic+1]].vN2,τ[ic+1].EQs[ic, ic, ess[ic,ic,ic+1]].vI2) + pi*Φ(τ[ic].EQs[ic, ic, ess[ic,ic,ic]].vN2,τ[ic].EQs[ic, ic, ess[ic,ic,ic]].vI2))+(1-P1)*((1-pi)*Φ(τ[ic+1].EQs[ic1, ic, ess[ic1,ic,ic+1]].vN2,τ[ic+1].EQs[ic1, ic, ess[ic1,ic,ic+1]].vI2)+pi*Φ(0,-K(c, cost)))))/(1-mp.β*(1-P1)*pi)
        vI2 = vN2 - K(c, cost)
        P2 = 0.

        τ[ic].EQs[ic1,ic,1] = EQ(P1, vN1, vI1, P2, vN2, vI2)
    end

    # if at (c,c2,c) edge
    for ic2 = 1:ic-1
        c2 = cost.C[ic2]
        vI2 = -K(c,cost) + mp.β*((1-pi)*Φ(τ[ic+1].EQs[ic, ic, ess[ic,ic,ic+1]].vN2,τ[ic+1].EQs[ic, ic, ess[ic,ic,ic+1]].vI2) + pi*Φ(τ[ic].EQs[ic, ic, ess[ic,ic,ic]].vN2,τ[ic].EQs[ic, ic, ess[ic,ic,ic]].vI2))
        
#        vN2search(x) = mp.β*((1 - pi)*Φ(τ[ic+1].EQs[ic, ic2,ess[ic,ic2,ic+1]].vN2, τ[ic+1].EQs[ic, ic2, ess[ic,ic2,ic+1]].vI2) + pi*Φ(x,vI2)) - x
#        vN2 = fzero(vN2search,0;order=2)

        #vN2 =  mp.β*vI2

        a = vI2
        A = mp.β*(1-pi)*Φ(τ[ic+1].EQs[ic,ic2,ess[ic,ic2,ic+1]].vN2, τ[ic+1].EQs[ic, ic2, ess[ic,ic2,ic+1]].vI2)
        B = mp.β*(1-pi)
        vN2 = 0
        if a == 0
            if  A >= 0
                vN2 = -A/(B-1)
            elseif A < 0
                vN2 = A
            end
        elseif a > 0
            if A >= a
                vN2 = -A/(B-1)
            elseif 0 < A < a
                if B <= (a-A)/a
                    vN2 = a*B+A
                elseif (a-A)/a < B
                    vN2 = -A/(B-1)
                end
            elseif A <= 0
                x = a*B+A
            end
        elseif a < 0
            if A >= 0
                vN2 = -A/(B-1)
            elseif a < A < 0
                if B <= (a-A)/a
                    vN2 = -A/(B-1)
                elseif (a-A)/a < B
                    vN2 = a*B+A
                end
            elseif A <= a
                vN2 = a*B+A
            end
        end

        P2 = 0.

        vN1 = (r1(c,c2) + mp.β*(P2*H1(ic,ic,ic, pi, τ, ess)+(1-P2)*((1-pi)*Φ(τ[ic+1].EQs[ic, ic2, ess[ic,ic2,ic+1]].vN1,τ[ic+1].EQs[ic, ic2, ess[ic,ic2,ic+1]].vI1)+pi*Φ(0,-K(c, cost)))))/(1-mp.β*(1-P2)*pi)
        vI1 = vN1-K(c, cost)
        P1 = 0.

        τ[ic].EQs[ic,ic2,1] = EQ(P1, vN1, vI1, P2, vN2, vI2)
    end
end

#=
function solve_interior!(cost::CostParameters, mp::RLSModel, ic::Int64, iC1::Int64, iC2::Int64, τ::Array{Stage,1}, ess::ESS)
    c = cost.C[ic]
    for ic2 in iC2:ic-1
        for ic1 in 1:ic-1
            find_interior!(mp, cost, ic1, ic2, ic, c, τ, ess)
        end
    end
end

solve_interior!(cost::CostParameters, mp::RLSModel, ic::Int64, τ, ess::ESS) = solve_interior!(cost, mp, ic, 1, 1, τ, ess)
 =#
function solve_interior!(cost::CostParameters, mp::RLSModel, ic::Int64, iC1::Int64, iC2::Int64, τ::Array{Stage,1}, ess::ESS)
            find_interior!(mp, cost, iC1, iC2, ic, cost.C[ic], τ, ess)
end
function solve_interior!(cost::CostParameters, mp::RLSModel, ic::Int64, τ::Array{Stage,1}, ess::ESS)
    for ic2 in 1:ic-1
        for ic1 in 1:ic-1
            find_interior!(mp, cost, ic1, ic2, ic, cost.C[ic], τ, ess)
        end
    end
end


function find_interior!(mp, cost, ic1, ic2, ic, c, τ::Array{Stage,1}, ess::ESS)
    pi = cost.π[ic]

    ## First player 1 -> leads to P2 candidates
    a = r1(cost.C[ic1], cost.C[ic2]) - K(c, cost) +mp.β*H1(ic, ic2, ic, pi, τ, ess)
    b = mp.β*(H1(ic, ic, ic, pi, τ, ess)-H1(ic, ic2, ic, pi, τ, ess))
    d = r1(cost.C[ic1],cost.C[ic2]) + mp.β*(1-pi)*Φ(τ[ic+1].EQs[ic1, ic2, ess[ic1,ic2,ic+1]].vN1,τ[ic+1].EQs[ic1, ic2, ess[ic1,ic2,ic+1]].vI1)
    e = mp.β*H1(ic1, ic, ic, pi, τ, ess)        - mp.β*(1-pi)*Φ(τ[ic+1].EQs[ic1, ic2, ess[ic1,ic2,ic+1]].vN1,τ[ic+1].EQs[ic1, ic2, ess[ic1,ic2,ic+1]].vI1)

    pa = - mp.β * pi * b
    pb = e + (mp.β*pi-1) * b - mp.β*pi*a
    pc = d + (mp.β*pi-1) * a

    pstar2 = quad(pa, pb, pc)

    ## Then player 2 -> leads to P1 candidates
    A = r2(cost.C[ic1], cost.C[ic2]) - K(c, cost) + mp.β*H2(ic1, ic, ic, pi, τ, ess)
    B = mp.β*(H2(ic, ic, ic, pi, τ, ess)-H2(ic1, ic, ic, pi, τ, ess))
    D = r2(cost.C[ic1],cost.C[ic2]) + mp.β*(1-pi)*Φ(τ[ic+1].EQs[ic1, ic2, ess[ic1,ic2,ic+1]].vN2,τ[ic+1].EQs[ic1, ic2, ess[ic1,ic2,ic+1]].vI2)
    E = mp.β*H2(ic, ic2, ic, pi, τ, ess)- mp.β*(1-pi)*Φ(τ[ic+1].EQs[ic1, ic2, ess[ic1,ic2,ic+1]].vN2,τ[ic+1].EQs[ic1, ic2, ess[ic1,ic2,ic+1]].vI2)

    qa = - mp.β * pi * B
    qb = E + (mp.β*pi-1) * B - mp.β*pi*A
    qc = D + (mp.β*pi-1) * A

    pstar1 = quad(qa, qb, qc)

    count = 0
    for i = 1:length(pstar1)
        for j = 1:length(pstar2)
            if i in 1:2 && j in 1:2
                exP1 = pc + pb * pstar2[j] + pa * pstar2[j]^2 < 0. ? 1. : 0.
                exP2 = qc + qb * pstar1[i] + qa * pstar1[i]^2 < 0. ? 1. : 0.

                if abs(exP1 - pstar1[i]) < 1e-7 && abs(exP2-pstar2[j]) < 1e-7
                    count += 1
                    vI1 = a + b*pstar2[j]
                    vN1 = (d + e*pstar2[j] + mp.β*pi*(1-pstar2[j])*(a+b*pstar2[j]))*pstar1[i]+(1-pstar1[i])*(d+e*pstar2[j])/(1-mp.β*pi*(1-pstar2[j]))
                    vI2 = A + B*pstar1[i]
                    vN2 = (D + E*pstar1[i] + mp.β*pi*(1-pstar1[i])*(A+B*pstar1[i]))*pstar2[j]+(1-pstar2[j])*(D+E*pstar1[i])/(1-mp.β*pi*(1-pstar1[i]))
                    τ[ic].EQs[ic1, ic2, count] = EQ(pstar1[i],vN1,vI1,pstar2[j],vN2,vI2)
                end
            elseif i > 2 && j > 2 && pstar1[i] >= 0 && pstar2[j] >= 0 && pstar1[i] <= 1 && pstar2[j] <= 1
                if pstar1[i] < 0 || pstar2[j] < 0 || pstar2[j] >1 || pstar1[i] > 1
                    println(pstar2[j])
                    println(pstar1[i])
                end
                count += 1
                v1 = a + b * pstar2[j]
                v2 = A + B * pstar1[i]
                @inbounds τ[ic].EQs[ic1, ic2, count] = EQ(pstar1[i],v1,v1,pstar2[j],v2,v2)
            end
        end
    end
#    #τ[ic].ne[ic1, ic2] = count
    ess.bases[1 + div(ess.nC*(ess.nC+1)*(2*ess.nC+1),6) - div(ic*(ic+1)*(2*ic+1),6) + 2*(ic - 1) + sub2ind((ic-1,ic-1),ic1,ic2)]=count
end