using Test 
using JuMP
using GLPK

default_optimizer = GLPK.Optimizer

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
        for i = 1:2  
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
    m = Model(default_optimizer)
    @testset "testing_max_of_reals" begin
        for i=1:100
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
    end
end

test_abs()
test_unit_step()
test_unit_step_times_real_var()
test_max_real()