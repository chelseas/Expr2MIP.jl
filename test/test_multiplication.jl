# testing
include("../src/parsing.jl")
include("../src/encodings.jl")
using Test
using Gurobi
using SymEngine
 using LazySets

# ENV["JULIA_DEBUG"] = Main
optimizer=Gurobi.Optimizer
function model()
    m = Model(optimizer)
    set_optimizer_attribute(m, "OutputFlag", 0)
    return m
end

# This script written to make sure the change of multiplication from OVERT  to bilinear work

function test1()
    # test optima from optimization vs analytical optima 
    # fun is a two argument, two variable function 
    fun(x,y) = x*y
    expr = Expr(fun(Basic("x"), Basic("y")))
    domain = Hyperrectangle(low=[-1., -1.], high=[1.,1.])
    
    m = model()
    x = @variable(m, base_name="x", lower_bound=low(domain)[1], upper_bound=high(domain)[1])
    y = @variable(m, base_name="y", lower_bound=low(domain)[2], upper_bound=high(domain)[2])
    params= EncodingParameters("interval", 1e-2, 1)
    d = Dict()
    con_ref, z = add_constraint!(m, expr, :z; params=params, expr_map=d)
    @objective(m, Max, z)
    optimize!(m)
    z_max_opt = value(z)
    # analytical max is 1 
    @test z_max_opt ≈ 1

    # now do min 
    @objective(m, Min, z)
    optimize!(m)
    z_min_opt = value(z)
    # analytical min is -1 
    @test z_min_opt ≈ -1.
end

function test2()
    # now more complex with constants and transcendentals
    fun(x,y) = sin(10*x*y) 
    expr = Expr(fun(Basic("x"), Basic("y")))
    domain = Hyperrectangle(low=[-1., -1.], high=[1.,1.])
    
    m = model()
    x = @variable(m, base_name="x", lower_bound=low(domain)[1], upper_bound=high(domain)[1])
    y = @variable(m, base_name="y", lower_bound=low(domain)[2], upper_bound=high(domain)[2])
    N=-1
    params= EncodingParameters("interval", 1e-2, N)
    d = Dict()
    con_ref, z = add_constraint!(m, expr, :z; params=params, expr_map=d)
    @objective(m, Max, z)
    optimize!(m)
    z_max_opt = value(z)
    # analytical max is 1 
    @test z_max_opt >= 1 # may be larger due to overapprox
    println("analytical opt is: 1, optimization opt with N=$N is: ", z_max_opt)

    # now do min 
    @objective(m, Min, z)
    optimize!(m)
    z_min_opt = value(z)
    # analytical min is -1 
    @test z_min_opt <= -1.
    println("analytical opt is: -1, optimization opt with N=$N is: ", z_min_opt)
end

test1()
test2()
