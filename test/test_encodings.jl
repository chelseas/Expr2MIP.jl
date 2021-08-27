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

function test_unit_step_times_real_var()
    m = Model(default_optimizer)
    @testset "testing_unit_step_times_var" begin
        for i = 1:100
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
            optimize!(m)
            @test value(z) ≈ (const_ẑ > 0 ? const_x : 0)
        end
    end
end

test_abs()
test_unit_step_times_real_var()