include("../src/parsing.jl")
using Test 
using Gurobi
using JuMP
using IntervalArithmetic

m = Model(Gurobi.Optimizer)
@test find_bounds(m, 0.0) == (0.0, 0.0)

# test interval arithmetic bounds 
function test_int_arith_bounds()
    model = Model(Gurobi.Optimizer)
    v = @variable(model, v[1:3], lower_bound=-1, upper_bound=1.0)
    expr = 1*v[1] + 2*v[2] + 3*v[3] - 12

    int_time = @elapsed find_bounds_through_int_arithmetic(model, expr) == interval(-18,-6)
    opt_time = @elapsed find_bounds_through_opt(model, expr) == interval(-18, -6)
    println("int time is $int_time and opt time is $opt_time")
end

function test_int_arith_bounds2()
    model = Model(Gurobi.Optimizer)
    v = @variable(model, v[1:3], lower_bound=-1, upper_bound=1.0)
    c = rand(4)
    expr = c[1]*v[1] + c[2]*v[2] + c[3]*v[3] - c[4]
    answer_min = sum(-1 .* c)
    answer_max = sum(c[1:3]) - c[4]
    answer = (answer_min, answer_max)

    int_time = @elapsed find_bounds_through_int_arithmetic(model, expr) 
    opt_time = @elapsed find_bounds_through_opt(model, expr)
    @assert all(find_bounds_through_int_arithmetic(model, expr) .≈ answer) 
    @assert all(find_bounds_through_opt(model, expr) .≈ answer)
    println("int time is $int_time and opt time is $opt_time")
    println("int time is: $(opt_time / int_time) times faster")
    return int_time, opt_time
end
