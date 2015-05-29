Base.run(`clear`)
println("Starting")
println("----------------------------")
using Leapfrog
println("JIT compilation")
mp, cost = init_all(1);
stags, ess, number_of_equilibria = solve(cost, mp);
println("----------------------------")
println("Checking # of equilibria")
println(number_of_equilibria)
println("----------------------------")
println("Timing code")
times = [@elapsed solve(cost,mp) for i = 1:10];
times_mean = mean(times)
times_min = minimum(times)
times_max = maximum(times)
times_var = std(times)
println("----------------------------")
@printf("Mean    Min    Max    Std\n")
@printf("%4.4f  %4.4f %4.4f %4.4f\n", times_mean, times_min, times_max, times_var)
println("----------------------------")
