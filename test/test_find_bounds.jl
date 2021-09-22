include("../src/parsing.jl")
using Test 
using Gurobi
using JuMP

m = Model(Gurobi.Optimizer)
@test find_bounds(m, 0.0) == (0.0, 0.0)
