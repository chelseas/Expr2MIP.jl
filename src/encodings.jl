using JuMP
global MAX_COUNT=0
global ABS_COUNT=0

function encode_abs!(model, input::VariableRef, LB, UB)
    """
    This function encodes absolute value as a mixed-integer program. LB and UB are lower and upper bounds on the input.
    """
    @assert LB <= UB
    output_var = @variable(model, lower_bound=0.0, base_name="t_$(ABS_COUNT)")
    δ = @variable(model, binary=true, base_name="δ_$(ABS_COUNT)")
    x⁺= @variable(model, lower_bound=0.0, base_name="x⁺_$(ABS_COUNT)")
    x⁻= @variable(model, lower_bound=0.0, base_name="x⁻_$(ABS_COUNT)")

    @constraint(model, input == x⁺ - x⁻)
    @constraint(model, output_var == x⁺ + x⁻)
    @constraint(model, x⁺ <= UB*δ)
    @constraint(model, x⁻ <= abs(LB)*(1-δ))
    global ABS_COUNT += 1

    return output_var
end

function encode_max_real!(model, inputs::Array{VariableRef}, LBs, UBs)
    # TODO
end

function encode_unit_step!(model, input::VariableRef, LB, UB)
    # TODO
end

function encode_unit_step_times_var!(model, ẑ::VariableRef, x::VariableRef, δ::VariableRef, l, u, γ, ζ)
    """
    This function encodes a unitstep*realvar using upper and lower bounds. l,u bound ẑ and γ, ζ bound x.
        l, u, γ, ζ are constants. 
    """
    # z = unit_step(ẑ)*x
    @assert JuMP.is_binary(δ)
    @assert l <= u 
    @assert γ <= ζ
    z = @variable(model, base_name="z") # output    
    @constraint(model, ẑ <= u*δ)
    @constraint(model, ẑ >= l*(1-δ))
    @constraint(model, z <= x - γ*(1-δ))
    @constraint(model, z >= x - ζ*(1-δ))
    @constraint(model, z >= γ*δ)
    @constraint(model, z <= ζ*δ)
    return z
end