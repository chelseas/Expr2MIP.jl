# test memoization
# pass a dictionary mapping expressions to JuMP stuff (Affine expressions or VariableRefs)

using Gurobi
using JuMP
include("../Expr2MIP.jl/src/parsing.jl")

function test1()
    model = Model(Gurobi.Optimizer)
    expr_map = Dict()
    x = @variable(model, base_name="x")
    set_lower_bound(x, 0.0)
    set_upper_bound(x, 1.0)
    add_constraint!(model, :(x), :e_1, expr_map=expr_map)
    print(model)
    add_constraint!(model, :(x), :e_2, expr_map=expr_map)
    print(model)
end

function test2()
    model = Model(Gurobi.Optimizer)
    expr_map = Dict()
    x = @variable(model, base_name="x")
    set_lower_bound(x, 0.0)
    set_upper_bound(x, 1.0)
    add_constraint!(model, :(2*x - 6), :e_1, expr_map=expr_map)
    print(model)
    add_constraint!(model, :(2*x - 6), :e_2, expr_map=expr_map)
    print(model)
end 

function test3()
    model = Model(Gurobi.Optimizer)
    expr_map = Dict()
    params = EncodingParameters("opt", 2e-2)
    x = @variable(model, base_name="x")
    set_lower_bound(x, 0.0)
    set_upper_bound(x, 1.0)
    add_constraint!(model, :(x*x), :e_1, expr_map=expr_map)
    model
    add_constraint!(model, :(x*x), :e_2, expr_map=expr_map)
    model

    add_constraint!(model, :(cos(x)), :e_3, expr_map=expr_map)
    model

    add_constraint!(model, :((x*x) + cos(x)), :e_4, expr_map=expr_map)
    model

    add_constraint!(model, :((x*x) * cos(x)), :e_5, expr_map=expr_map)
    model
end

