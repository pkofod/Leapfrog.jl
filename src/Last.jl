function solve_last_corner!(cost::CostParameters, mp::RLSModel, τ)
    c = cost.Cmin
    P1 = 0.
    vN1 = (mp.β*Φ(0, -K(c, cost)))/(1 - mp.β)
    vI1 = vN1 - K(c, cost)
    vN2 = (mp.β*Φ(0, -K(c, cost)))/(1 - mp.β)
    vI2 = vN2 - K(c, cost)
    P2 = 0.
    τ[end].EQs[end, end, 1] = EQ(P1, vN1, vI1, P2, vN2, vI2)
#    τ[end].ne[end, end] = 1
end

function solve_last_edge!(cost::CostParameters, mp::RLSModel, τ)
    nC = cost.nC
    c = cost.Cmin

    # if at (c1,c,c) edge
    for ic1 = 1:nC-1
        c1 = cost.C[ic1]
        vI1 =  - K(0., cost) + mp.β*Φ(τ[end].EQs[end, end, 1].vN1,τ[end].EQs[end, end, 1].vI1)
        #vN1search(x) =  mp.β*Φ(x,vI1)-x
        #vN1 = fzero(vN1search,0;order=2)
        #vN1 =  mp.β*vI1
        vN1 = 0. # Since vI1 in the corner is negative!
        P1 = 0.
        vN2 = (r2(c1, c) + mp.β*(P1*0+(1-P1)*Φ(0,-K(c, cost))))/(1-mp.β*(1-P1))
        vI2 = vN2 - K(c, cost)
        P2 = 0.

        τ[end].EQs[ic1, end, 1] = EQ(P1, vN1, vI1, P2, vN2, vI2)
#        τ[end].ne[ic1, end] = 1
    end
    for ic2 = 1:nC-1
        c2 = cost.C[ic2]
        vI2 = r2(c, c2) - K(c, cost) + mp.β*Φ(τ[end].EQs[end, end, 1].vN2,τ[end].EQs[end, end, 1].vI2)
        #vN2search(x) = r2(c, c2) + mp.β*Φ(x,vI2)-x
        #vN2 = fzero(vN2search,0;order=2)
        #vN2 =  mp.β*vI2
        vN2 = 0.
        P2 = 0.
        vN1 = (r1(c, c2) + mp.β*(P2*0+(1-P2)*Φ(0, -K(c, cost))))/(1 - mp.β*(1-P2))
        vI1 = vN1-K(c, cost)
        P1 = 0.

        τ[end].EQs[end, ic2, 1] = EQ(P1, vN1, vI1, P2, vN2, vI2)
#        τ[end].ne[end, ic2] = 1
    end
end


function solve_last_interior!(cost::CostParameters, mp::RLSModel, ic1::Int64, ic2::Int64, τ, ess)
            find_interior_last!(cost, mp, ic1, ic2, ess.nC, cost.C[ess.nC], τ, ess)
end

function solve_last_interior!(cost::CostParameters, mp::RLSModel, τ, ess)
    for ic2 in 1:ess.nC-1
        for ic1 in 1:ess.nC-1
            find_interior_last!(cost, mp, ic1, ic2, ess.nC, cost.C[ess.nC], τ, ess)
        end
    end
end
#=
function solve_last_interior!(cost::CostParameters, mp::RLSModel, τ, ESS)
    nC = cost.nC
    ic = nC
    c = cost.C[end]
    for ic2 in 1:nC-1
        for ic1 in 1:nC-1
            find_interior_last!(cost, mp, ic1, ic2, ic, c, τ, ESS)
        end
    end
end
=#
function find_interior_last!(cost::CostParameters, mp::RLSModel, ic1, ic2, ic, c, τ, ESS)
    pi = cost.π[ic]

    ## First player 1 -> leads to P2 candidates
    a = r1(cost.C[ic1], cost.C[ic2]) - K(c, cost) + mp.β*H1(ic, ic2, pi, τ, ESS) # OK
    b = mp.β*(H1(ic, ic, pi, τ, ESS)-H1(ic, ic2, pi, τ, ESS)) # OK
    d = r1(cost.C[ic1],cost.C[ic2])
    e = mp.β*H1(ic1, ic, pi, τ, ESS)

    pa = - mp.β * pi * b
    pb = mp.β*H1(ic1, ic, pi, τ, ESS) + (mp.β*pi-1)*b - mp.β*pi*a
    pc = r1(cost.C[ic1],cost.C[ic2]) + (mp.β*pi-1) * a

    pstar2 = quad(pa, pb, pc)

    ## Then player 2 -> leads to P1 candidates
    A = r2(cost.C[ic1], cost.C[ic2]) - K(c, cost) + mp.β*H2(ic1, ic, pi, τ, ESS) #OK
    B = mp.β*(H2(ic, ic, pi, τ, ESS)-H2(ic1, ic, pi, τ, ESS))#OK
    D = r2(cost.C[ic1],cost.C[ic2])
    E = mp.β*H2(ic, ic2, pi, τ, ESS)

    qa = - mp.β * pi * B
    qb = mp.β*H2(ic, ic2, pi, τ, ESS) + (mp.β*pi-1) * B - mp.β*pi*A
    qc = r2(cost.C[ic1],cost.C[ic2]) + (mp.β*pi-1) * A
    pstar1 = quad(qa, qb, qc)

    count = 0
    for i = 1:length(pstar1)
        for j = 1:length(pstar2)
            if i in 1:2 && j in 1:2
                exP1 = pc + pb * pstar2[j] + pa * pstar2[j]^2 < 0. ? 1. : 0.
                exP2 = qc + qb * pstar1[i] + qa * pstar1[i]^2 < 0. ? 1. : 0.

                if abs(exP1 - pstar1[i]) < 1e-8 && abs(exP2-pstar2[j]) < 1e-8
                    count += 1
                    vI1 = a + b*pstar2[j]
                    vN1 = (d + e*pstar2[j] + mp.β*pi*(1-pstar2[j])*(a+b*pstar2[j]))*pstar1[i]+(1-pstar1[i])*(d+e*pstar2[j])/(1-mp.β*pi*(1-pstar2[j]))
                    vI2 = A + B*pstar1[i]
                    vN2 = (D + E*pstar1[i] + mp.β*pi*(1-pstar1[i])*(A+B*pstar1[i]))*pstar2[j]+(1-pstar2[j])*(D+E*pstar1[i])/(1-mp.β*pi*(1-pstar1[i]))

                    τ[end].EQs[ic1, ic2, count] = EQ(pstar1[i],
                                                     vN1,
                                                     vI1,
                                                     pstar2[j],
                                                     vN2,
                                                     vI2)
                end
            elseif i > 2 && j > 2 && pstar1[i] >= 0 && pstar2[j] >= 0 && pstar1[i] <= 1 && pstar2[j] <= 1
                count += 1
                v1 = a + b * pstar2[j]
                v2 = A + B * pstar1[i]
                @inbounds τ[end].EQs[ic1, ic2, count] = EQ(pstar1[i],v1,v1,pstar2[j],v2,v2)
            end
        end
    end
    ESS.bases[1 + div(ESS.nC*(ESS.nC+1)*(2*ESS.nC+1),6) - div(ic*(ic+1)*(2*ic+1),6) + 2*(ic - 1) + sub2ind((ic-1,ic-1),ic1,ic2)]=count
end