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

    return output_var::VariableRef
end

function get_u_max_i(us::Array, i)
    """Helper function for encode_max_real!
    # get max excluding ith value
    """
    return maximum(us[[1:(i-1)...,(i+1):end...]])
end

function encode_max_real!(model, inputs::Array{VariableRef}, LBs::Array, UBs::Array)
    """
    This function encodes the maximum of several real-valued variables. 
        Each input is assumed to have both a lower bound and upper bound, contained respectively in the arrays LBs and UBs
    """
    @assert all(LBs .<= UBs) #assert each pairing of (LB, UB) has LB <= UB
    l_max = maximum(LBs)
    # eliminate any x_i from consideration where u_i <= l_max since we know y >= l_max >= u_i >= x_i
    indices_to_keep = findall(u -> u > l_max, UBs)
    inputs = inputs[indices_to_keep]
    LBs = LBs[indices_to_keep]
    UBs = UBs[indices_to_keep]
    new_n = length(inputs)

    # create output variable
    y = @variable(model, base_name="y_$(MAX_COUNT)", lower_bound=l_max, upper_bound=maximum(UBs)) # y = max(x_1, x_2, ...)
    δeltas = @variable(model, [1:new_n], binary=true, base_name="δ_max_$(MAX_COUNT)")
    for i = 1:new_n
        u_max_i = get_u_max_i(UBs, i)
        xᵢ = inputs[i]
        aᵢ = δeltas[i]
        lᵢ = LBs[i]
        @constraint(model, y <= xᵢ + (1 - aᵢ)*(u_max_i - lᵢ))
        @constraint(model, y >= xᵢ)
    end
    @constraint(model, sum(δeltas) == 1)
    return y::VariableRef
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
    return z::VariableRef
end