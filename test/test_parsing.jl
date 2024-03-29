# testing
include("../src/parsing.jl")
include("../src/encodings.jl")
using Test
using Gurobi
# ENV["JULIA_DEBUG"] = Main
optimizer=Gurobi.Optimizer
model() = Model(with_optimizer(optimizer, OutputFlag=0))

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
    m = model()
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

        function f1()
            m = model()
            in_ref = @variable(m, x, lower_bound=-10.2, upper_bound=6.7)
            con_ref, ovar_ref = add_constraint!(m, :(abs(x)), :o)
            print(m)
            optimize!(m)
            @test value(ovar_ref) >= 0
            @constraint(m, x == -5)
            optimize!(m)
            @test abs(value(x)) == value(ovar_ref)
        end

        function f2()
            m = model()
            in_ref = @variable(m, x, lower_bound=-1.1, upper_bound=1.1)
            con_ref, ovar_ref = add_constraint!(m, :(x + abs(x)), :o)
            print(m)
            const_val = rand()*2-1
            @constraint(m, in_ref == const_val)
            optimize!(m)
            @test value(ovar_ref) == const_val + abs(const_val)
        end

        function f3()
            m = model()
            in_ref = @variable(m, x, lower_bound=-1.1, upper_bound=1.1)
            add_constraint!(m, :(x + 5*abs(x)), :o)
            print(m)
        end

        function f4()
            m = model()
            in_ref = @variable(m, x, lower_bound=-1.1, upper_bound=1.1)
            con_ref, ovar_ref = add_constraint!(m, :(5*x + 5*abs(x)), :o)
            print(m)
            const_val = rand()*2-1
            @constraint(m, in_ref == const_val)
            optimize!(m)
            @test value(ovar_ref) == 5*const_val + 5*abs(const_val)
        end

        function f5()
        m = model()
        x = @variable(m, x, lower_bound=-1.1, upper_bound=1.1)
        y = @variable(m, y, lower_bound=-1.1, upper_bound=1.1)
        add_constraint!(m, :(abs(x) + abs(y)), :o)
        print(m)
        end

        function f6()
        m = model()
        in_ref = @variable(m, x, lower_bound=-1.1, upper_bound=1.1)
        add_constraint!(m, :(abs(-5 + 10*abs(x))), :o)
        print(m)
        end

        function f7()
        m = model()
        x = @variable(m, x, lower_bound=-1.1, upper_bound=1.1)
        y = @variable(m, y, lower_bound=-1.1, upper_bound=1.1)
        add_constraint!(m, :(4*x + 5*y - 62), :o)
        print(m)
        end

        function f8()
            m = model()
            ẑ = @variable(m, ẑ, lower_bound=-6.7, upper_bound=5.6)
            con_ref, ovar_ref = add_constraint!(m, :(unit_step(ẑ)+1), :o)
            print(m)
            const_ẑ = rand()*2 - 1
            @constraint(m, ẑ == const_ẑ)
            optimize!(m)
            @test value(ovar_ref) == (sign(const_ẑ)+1)/2 + 1
        end

        function f9()
            m = model()
            ẑ = @variable(m, ẑ, lower_bound=-6.7, upper_bound=5.6)
            x = @variable(m, x, lower_bound=-12.3, upper_bound=0.5)
            params = EncodingParameters()
            con_ref, ovar_ref = add_constraint!(m, :(unit_step(ẑ)*x), :o, params=params)
            print(m)
            const_ẑ = rand()*12.3 - 6.7
            const_x = rand()*12.8 - 12.3
            @constraint(m, ẑ == const_ẑ)
            @constraint(m, x == const_x)
            optimize!(m)
            @test value(ovar_ref) ≈ unit_step(const_ẑ)*const_x
        end

        function f10()
            m = model()
            x_ref = @variable(m, x, lower_bound=-1, upper_bound=1.0)
            y_ref = @variable(m, y, lower_bound=-1., upper_bound=1.0)
            params = EncodingParameters()
            con_ref, ovar_ref = add_constraint!(m, :(max(x,y)), :o, params=params)
            print(m)
            const_x = rand()*2 - 1
            const_y = rand()*2 - 1
            @constraint(m, x_ref == const_x)
            @constraint(m, y_ref == const_y)
            optimize!(m)
            @test value(ovar_ref) == max(const_x, const_y)
        end

        function f11()
            m = model()
            x_ref = @variable(m, x, lower_bound=-1, upper_bound=1.0)
            y_ref = @variable(m, y, lower_bound=-1., upper_bound=1.0)
            con_ref, ovar_ref = add_constraint!(m, :(min(x,y)), :o)
            print(m)
            const_x = rand()*2 - 1
            const_y = rand()*2 - 1
            @constraint(m, x_ref == const_x)
            @constraint(m, y_ref == const_y)
            optimize!(m)
            @test value(ovar_ref) == min(const_x, const_y)
        end

    @testset "testing encoding" begin
        f1()
        f2()
        f3()
        f4()
        f5()
        f6()
        f7()
        f8()
        f9()
        f10()
        f11()
    end
end
testing_encoding()

function test_symbols()
    m = model()
    x_ref = @variable(m, x, lower_bound=-1, upper_bound=1.0)
    con_ref, ovar_ref = add_constraint!(m, :(pi*(1/180)*x), :o)
    print(m)
end

function test_overt()
    m = Model(with_optimizer(Gurobi.Optimizer, OutputFlag=0))
    x_ref = @variable(m, x, lower_bound=-1, upper_bound=1.0)
    e = :(sin(x/2))
    con_ref, ovar_ref = add_constraint!(m, e, :o)
    # then fix the value of x, max and min oa.output and assert that this is an interval capturing the true value
    model = m;
    @constraint(model, x_ref == pi/4)
    @objective(model, Max, ovar_ref)
    optimize!(model)
    e_max = value(ovar_ref)
    @objective(model, Min, ovar_ref)
    optimize!(model)
    e_min = value(ovar_ref)
    e_true = sin(value(x_ref)/2) 
    @assert (e_min <= e_true) && (e_true <= e_max)
    println("(e_min, e_true, e_max) = ($(e_min), $(e_true), $(e_max))")
end

function test_overt_complex()
    # recommend plotting in wolfram alpha to visualize
    m = Model(with_optimizer(Gurobi.Optimizer, OutputFlag=0))
    x_ref = @variable(m, x, lower_bound=-10, upper_bound=10.0)
    e = :(max(sin(x/2), x))
    con_ref, ovar_ref = add_constraint!(m, e, :o)
    # then fix the value of x, max and min oa.output and assert that this is an interval capturing the true value
    model = m;
    # @constraint(model, x_ref == pi/4)
    @objective(model, Max, ovar_ref)
    optimize!(model)
    e_max = value(ovar_ref)
    @objective(model, Min, ovar_ref)
    optimize!(model)
    e_min = value(ovar_ref)
    e_true = max(sin(value(x_ref)/2), value(x_ref)) 
    @assert (e_min <= e_true) && (e_true <= e_max)
    println("(e_min, e_true, e_max) = ($(e_min), $(e_true), $(e_max))")
end

function test_max_rewrite()
    # test the re-write of max to relu
    # max(a,b) = 0.5*(a + b + relu(a-b) + relu(b-a))
    relu(x) = max(x, 0)
    x = rand()
    y = rand()
    @assert max(x,y) == 0.5*(x + y + relu(x - y) + relu(y - x))
end