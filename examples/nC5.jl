Base.run(`clear`)
println("Starting")
println("----------------------------")
using Leapfrog
println("JIT compilation")
@elapsed 1
mp_burn, cost_burn = init_all(1)
solve(cost_burn, mp_burn)
mp, cost = init_all(2);
@elapsed stags, ess, number_of_equilibria = solve(cost, mp)
println("----------------------------")
println("Checking # of equilibria")
println(number_of_equilibria)
println("----------------------------")