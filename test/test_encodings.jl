using Test 
using JuMP
#using GLPK
using Gurobi
#ENV["JULIA_DEBUG"] = Main

include("../src/encodings.jl")

default_optimizer = Gurobi.Optimizer #GLPK.Optimizer

function test_abs()
    m = Model(default_optimizer)
    @testset "testing_abs" begin
        for i = 1:100
            input = @variable(m)
            output = encode_abs!(m, input, -5., 5.)
            const_val = rand()*9 - 4.5 # in range -4.5, 4.5
            @constraint(m, input == const_val)
            optimize!(m)
            @test value(output) ≈ abs(const_val)
        end
    end    
end

function test_unit_step()
    m = Model(default_optimizer)
    @testset "testing_unit_step" begin
        for i = 1:100
            input = @variable(m)
            output = encode_unit_step!(m, input, -5., 5.)
            const_val = rand()*9 - 4.5 # in range -4.5, 4.5
            @constraint(m, input == const_val)
            optimize!(m)
            @test value(output) == Int(const_val > 0)
        end
    end 
end

function test_unit_step_times_real_var()
    m = Model(default_optimizer)
    @testset "testing_unit_step_times_var" begin
        for i = 1:100
            m = Model(default_optimizer)          
            ẑ = @variable(m, base_name="ẑ")
            x = @variable(m, base_name="x")
            δ = @variable(m, base_name="δ", binary=true)
            l, u = [-5., 5.]
            γ, ζ =  [-6., 34]
            z = encode_unit_step_times_var!(m, ẑ, x, δ, l, u, γ, ζ)

            const_ẑ = rand()*9 - 4.5 # in rangeẑ = @variable(m, base_name="ẑ")
            const_x = rand()*39 - 5.5 # in range -5.5, 33.5
            @constraint(m, ẑ == const_ẑ)
            @constraint(m, x == const_x)
            @show num_constraints(m, VariableRef, MOI.GreaterThan{Float64})
            @show num_constraints(m, VariableRef, MOI.ZeroOne)
            @show num_constraints(m, AffExpr, MOI.LessThan{Float64})
            optimize!(m)
            @show termination_status(m)

            @test value(z) ≈ (const_ẑ > 0 ? const_x : 0)
        end
    end
end

function test_max_real()
    @testset "testing_max_of_reals" begin
        for i=1:100
            m = Model(default_optimizer)
            inputs = @variable(m, [1:10], base_name="inputs")
            LBs = rand([1:3...], 10) # 10 random ints btw 1 and 3
            UBs = rand([0:3...], 10) .+ LBs 
            @assert all(LBs .<= UBs) #assert each pairing of (LB, UB) has LB <= UB
            y = encode_max_real!(m, inputs, LBs, UBs)
            # now fix variables and assert that the max is correct
            max_xᵢ = -Inf
            for i = 1:10
                xᵢ = rand()*(UBs[i] - LBs[i]) + LBs[i]
                @assert LBs[i] <= xᵢ <= UBs[i]
                max_xᵢ = max(max_xᵢ, xᵢ)
                @constraint(m, inputs[i] == xᵢ)
            end
            optimize!(m)
            @test value(y) == max_xᵢ
        end

        # edge case 
        m = Model(default_optimizer)
        x1 = @variable(m, base_name="x1")
        x2 = @variable(m, base_name="x2")
        x3 = @variable(m, base_name="x3")
        LBs = [0.0, -2.999999, -4.]
        UBs = [-0.0, 1.0, 0.9]
        y = encode_max_real!(m, [x1, x2, x3], LBs, UBs)
        @constraint(m, x2 == -1.0)
        @constraint(m, x3 == -.5)
        optimize!(m)
        @test value(y) == 0.0

    end
end

function test_relu()
    @testset "test_relus" begin
        for i = 1:100
            m = Model()
            set_optimizer(m, optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
            lb = randn()
            ub = lb + rand()
            println("bounds are: [$lb, $ub]")
            ẑ = @variable(m, base_name="relu input", lower_bound=lb, upper_bound=ub)
            print(m)
            z = encode_relu!(m, ẑ::t, lb, ub; relax="none")
            print(m)
            optimize!(m)
            relu(x) = max(x,0.0)
            @test value(z) == relu(value(ẑ))
        end
    end
end

function test_relu_rewrite_max()

    # NOTE: this will fail if the encode_relu! parameter relax defaults to "triangle" and that is expected. 
    @testset "test relu rewrite" begin
        for i=1:100
            LBs = randn(2)
            UBs = rand(2) .+ LBs 
            @assert all(LBs .<= UBs) #assert each pairing of (LB, UB) has LB <= UB

            m = Model(default_optimizer)
            x1 = @variable(m, base_name="x1", lower_bound=LBs[1], upper_bound=UBs[1])
            x2 = @variable(m, base_name="x2", lower_bound=LBs[2], upper_bound=UBs[2])
            args = [:(x1), :(x2)]
            print(m)
            m1 = encode!(m, Sym_f{:max}(), args::Array; bound_type="interval", relu_rewrite=false)
            print(m)
            # TODO: have the following version be smarter and not add a bunch of constraints if a max isn't needed.
            m2 = encode!(m, Sym_f{:max}(), args::Array; bound_type="interval", relu_rewrite=true)
            print(m)
            # now fix variables and assert that the max is correct
            max_xᵢ = -Inf
            for i = 1:2
                xᵢ = rand()*(UBs[i] - LBs[i]) + LBs[i]
                @assert LBs[i] <= xᵢ <= UBs[i]
                max_xᵢ = max(max_xᵢ, xᵢ)
                @constraint(m, JuMP.variable_by_name(m, "x$i") == xᵢ)
            end
            optimize!(m)
            @test value(m1) ≈ max_xᵢ
            @test value(m2) ≈ max_xᵢ
        end 
    end
end

test_abs()
test_unit_step()
test_unit_step_times_real_var()
test_max_real()