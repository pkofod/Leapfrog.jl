# Leapfrog.jl

This package solves the dynamic Bertrand model with cost-reducing investments in Ishkakov, Rust, and Schjerning (2015): "The Dynamics of Bertrand Competition with Cost-reducing Investments". Currently, the simultaneous move version of the game is implemented.

# Use

To solve a model, you have to load the package, create a variable of type CostParameter, and a variable of type RLSModel. Then you just have to call `solve`.

```julia
using Leapfrog
cost, mp = init_all(1) # 1 for the nC = 4 version
stages, ess, number_of_equilibria = solve(cost, mp)
```

A lot of this could be cleaner, but it works. Parameters of the model can be changed in Leapfrog.jl.
