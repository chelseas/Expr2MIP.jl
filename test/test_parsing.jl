# testing
include("../src/parsing.jl")
using Test

function test_sym_f()
    function test(a::Sym_f{:abs})
        println("a = $a")
    end
    function test(b::Sym_f{:exp})
        println("b = $b")
    end
    test(Sym_f(:abs))
    test(Sym_f(:exp))
end
test_sym_f()

function affine_test()
    @testset "test_affine" begin 
        @test is_affine(:(x+1))
        @test is_affine(:(-x+(2y-z)))
        @test !is_affine(:(log(x)))  
        @test !is_affine(:(x + x*z))  
        @test is_affine(:(x/6))      
        @test !is_affine(:(6/x))     
        @test is_affine(:(5*x))      
        @test is_affine(:(log(2)*x)) 
        @test is_affine(:(-x))       
    end
end
affine_test()

using Gurobi
function test_convert_affine_to_jump()
    m = JuMP.Model(Gurobi.Optimizer)
    e1 = :(3*(x + y) - z)
    @variable(m, base_name="x", lower_bound=0., upper_bound=1.3)
    @variable(m, base_name="y", lower_bound=0., upper_bound=1.3)
    @variable(m, base_name="z", lower_bound=0., upper_bound=1.3)
    jump_e1 = convert_affine_to_jump(e1, m)
    println(jump_e1)
    # @constraint(m, e1 == 4) -> does not  work
    @constraint(m, jump_e1 == 4) # works
    print(m)
end
test_convert_affine_to_jump()

function test_convert_step_times_var()
    @testset "test convert step times var" begin
    #println(convert_step_times_var(:(unit_step(z)*x)))
    @test convert_step_times_var(:(unit_step(z)*x)) == :(unit_step_times_var(z, x))
    #println(convert_step_times_var(:(5*unit_step(z)*x)))
    @test convert_step_times_var(:(5*unit_step(z)*x)) == :(5*unit_step_times_var(z, x))
    #println(convert_step_times_var(:5))
    @test convert_step_times_var(:5) == 5
    #println(convert_step_times_var(:5*6*7))
    @test convert_step_times_var(:5*6*7) == 5*6*7
    #println(convert_step_times_var(:(5x + 3y - 10)))
    @test convert_step_times_var(:(5x + 3y - 10)) == :(5x + 3y - 10)
    #println(convert_step_times_var(:(unit_step(5x + 3y - 10))))
    @test convert_step_times_var(:(unit_step(5x + 3y - 10))) == :(unit_step(5x + 3y - 10))
    #println(convert_step_times_var(:(unit_step(5x + 3y - 10)*q*5)))
    @test convert_step_times_var(:(unit_step(5x + 3y - 10)*q*5)) == :(5*unit_step_times_var(5x + 3y - 10, q))
    #println(convert_step_times_var(:(unit_step(5x + 3y - 10)*q*5 - 10*x)))
    @test convert_step_times_var(:(unit_step(5x + 3y - 10)*q*5 - 10*x)) == :(5*unit_step_times_var(5x + 3y - 10, q) - 10*x)
    #println(:(unit_step(unit_step(5x)*q*5 - 10*x)*z))
    #println(convert_step_times_var(:(unit_step(unit_step(5x)*q*5 - 10*x)*z)))
    @test convert_step_times_var(:(unit_step(unit_step(5x)*q*5 - 10*x)*z)) == :(unit_step_times_var(5*unit_step_times_var(5x, q) - 10*x, z))
    end
end
test_convert_step_times_var()

# testing encoding functions
function testing_encoding()
    @testset "testing encoding" begin
        m = Model(Gurobi.Optimizer)
        in_ref = @variable(m, x)
        con_ref, ovar_ref = add_constraint!(m, :(abs(x)), :o)
        print(m)
        optimize!(m)
        @test value(ovar_ref) >= 0
        @constraint(m, x == -5)
        optimize!(m)
        @test abs(value(x)) == value(ovar_ref)

        m = Model(Gurobi.Optimizer)
        in_ref = @variable(m, x)
        con_ref, ovar_ref = add_constraint!(m, :(x + abs(x)), :o)
        print(m)
        const_val = rand()*2-1
        @constraint(m, in_ref == const_val)
        optimize!(m)
        @test value(ovar_ref) == const_val + abs(const_val)

        m = Model(Gurobi.Optimizer)
        add_constraint!(m, :(x + 5*abs(x)), :o)
        print(m)

        m = Model(Gurobi.Optimizer)
        add_constraint!(m, :(5*x + 5*abs(x)), :o)
        print(m)

        m = Model(Gurobi.Optimizer)
        add_constraint!(m, :(abs(x) + abs(y)), :o)
        print(m)

        m = Model(Gurobi.Optimizer)
        add_constraint!(m, :(abs(-5 + 10*abs(x))), :o)
        print(m)

        m = Model(Gurobi.Optimizer)
        add_constraint!(m, :(4*x + 5*y - 62), :o)
        print(m)

        m = Model(Gurobi.Optimizer)
        ẑ = @variable(m, ẑ)
        con_ref, ovar_ref = add_constraint!(m, :(unit_step(ẑ)), :o)
        print(m)
        const_ẑ = rand()*2 - 1
        @constraint(m, ẑ == const_ẑ)
        optimize!(m)
        @test value(ovar_ref) == (sign(const_ẑ)+1)/2

        m = Model(Gurobi.Optimizer)
        add_constraint!(m, :(unit_step(z)*x), :o)
        print(m)

        m = Model(Gurobi.Optimizer)
        x_ref = @variable(m, x)
        y_ref = @variable(m, y)
        con_ref, ovar_ref = add_constraint!(m, :(max(x,y)), :o)
        print(m)
        const_x = rand()*2 - 1
        const_y = rand()*2 - 1
        @constraint(m, x_ref == const_x)
        @constraint(m, y_ref == const_y)
        optimize!(m)
        @test value(ovar_ref) == max(const_x, const_y)

        m = Model(Gurobi.Optimizer)
        x_ref = @variable(m, x)
        y_ref = @variable(m, y)
        con_ref, ovar_ref = add_constraint!(m, :(min(x,y)), :o)
        print(m)
        const_x = rand()*2 - 1
        const_y = rand()*2 - 1
        @constraint(m, x_ref == const_x)
        @constraint(m, y_ref == const_y)
        optimize!(m)
        @test value(ovar_ref) == min(const_x, const_y)
    end
end
testing_encoding()
