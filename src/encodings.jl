using JuMP
global MAX_COUNT=0
global ABS_COUNT=0
global RELU_COUNT=0
global UNIT_STEP_COUNT=0
global M = 1000000 # very large number.
global NEW_OVERT_VAR_COUNT=0

using NeuralVerification
using NeuralVerification: encode_relu, BoundedMixedIntegerLP, TriangularRelaxedLP

# These functions take a JuMP-compatible type, which is either a VariableRef or 
# GenericAffExpr and return a VariableRef corresponding to the output variable

t = Union{VariableRef, GenericAffExpr, AffExpr, Real}

function encode_abs!(model, input::t, LB, UB)
    """
    This function encodes absolute value as a mixed-integer program. LB and UB are lower and upper bounds on the input.
    """
    @assert LB <= UB
    output_var = @variable(model, lower_bound=0.0, upper_bound=max(-LB, UB), base_name="t_$(ABS_COUNT)")
    δ = @variable(model, binary=true, base_name="δ_abs_$(ABS_COUNT)")
    x⁺= @variable(model, lower_bound=0.0, upper_bound=max(UB, 0.0), base_name="x⁺_$(ABS_COUNT)")
    x⁻= @variable(model, lower_bound=0.0, upper_bound=max(-LB, 0.0), base_name="x⁻_$(ABS_COUNT)")

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

function check_if_need_max(LBs, UBs)
    """
    A helper function to determine if we actually need to take a max or not.
    """
    @debug("LBs: $(LBs)")
    @debug("UBs: $(UBs)")
    @assert all(LBs .<= UBs) #assert each pairing of (LB, UB) has LB <= UB
    l_max = maximum(LBs)
    # eliminate any x_i from consideration where u_i < l_max since we know y >= l_max >= u_i >= x_i
    indices_to_keep = findall(u -> u >= l_max, UBs)
    need_max = length(indices_to_keep) > 1
    return need_max, indices_to_keep, l_max
end

function encode_max_real!(model, inputs::Array, LBs::Array, UBs::Array)
    """
    This function encodes the maximum of several real-valued variables. 
        Each input is assumed to have both a lower bound and upper bound, contained respectively in the arrays LBs and UBs
    """
    need_max, indices_to_keep, l_max = check_if_need_max(LBs, UBs)
    u_max = maximum(UBs)
    inputs = inputs[indices_to_keep]
    LBs = LBs[indices_to_keep]
    UBs = UBs[indices_to_keep]
    new_n = length(inputs)
    @assert new_n >= 1
    #println("new_n = $new_n")
    @debug("Performing a max over $inputs")

    if !need_max # (new_n == 1)
        # there is no need to actually take a max 
        y = inputs[1] # should just be 1, the variable that is alwa@debug "there is no need to actually take a max"ys the max
        @debug "there is no need to actually take a max." y UBs
    else # need max
        # create output variable
        y = @variable(model, base_name="y_$(MAX_COUNT)", lower_bound=l_max, upper_bound=u_max) # y = max(x_1, x_2, ...)
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
        global MAX_COUNT += 1
    end
    return y::t
end

function encode_unit_step!(model, input::t, LB, UB)
    """
    This function encodes a unit step all by itself. e.g. δ = unit_step(x)
    """
    @assert LB <= UB 
    δ = @variable(model, binary=true, base_name="unit_step_$(UNIT_STEP_COUNT)", lower_bound=0.0, upper_bound=1.0)
    global UNIT_STEP_COUNT += 1
    @constraint(model, UB*δ >= input)
    @constraint(model, input >= LB*(1 - δ))
    return δ
end

function encode_unit_step_times_var!(model, ẑ::t, x::t, δ::VariableRef, l, u, γ, ζ)
    """
    This function encodes a unitstep*realvar using upper and lower bounds. l,u bound ẑ and γ, ζ bound x.
        l, u, γ, ζ are constants. 
    """
    # z = unit_step(ẑ)*x
    @assert JuMP.is_binary(δ)
    @assert l <= u 
    @assert γ <= ζ
    always_off = u < 0
    if always_off
        # output variable 
        @debug "Always off"
        z = @variable(model, base_name="z", lower_bound=0.0, upper_bound=0.0)            
    else          
        # output variable 
        z = @variable(model, base_name="z", lower_bound=min(0.0, γ), upper_bound=max(0.0, ζ))
        @constraint(model, ẑ <= u*δ)
        @constraint(model, ẑ >= l*(1-δ))
        @constraint(model, z <= x - γ*(1-δ))
        @constraint(model, z >= x - ζ*(1-δ))
        @constraint(model, z >= γ*δ)
        @constraint(model, z <= ζ*δ)
    end
    return z::VariableRef
end

function encode_relu!(model, ẑ::t, lower, upper; relax="none")
    """
    This function encodes the z = relu(ẑ) aka max(ẑ,0).
        relax can be:
        - "none"
        - "triangle"
        - ? 
    """
    # thoughts: could I call a zonotope method here?
    # start by adding some cases for different kinds of encoding. exact and different relaxations.  
    println("Encoding relu using $(relax) relaxation. lower bound: $lower, upper bound $upper")
    if relax == "none"
        z = @variable(model, base_name="z_relu_$(RELU_COUNT)", lower_bound=max(0.0, lower), upper_bound=max(0.0, upper))
        δ = @variable(model, base_name="δ_relu_$(RELU_COUNT)", binary=true)
        global RELU_COUNT += 1
        NeuralVerification.encode_relu(BoundedMixedIntegerLP(), model, ẑ, z, δ, lower, upper)
        return z::t
    elseif relax == "triangle"
        z = @variable(model, base_name="z_relu_$(RELU_COUNT)", lower_bound=max(0.0, lower), upper_bound=max(0.0, upper))
        global RELU_COUNT += 1
        NeuralVerification.encode_relu(TriangularRelaxedLP(), model, ẑ, z, lower, upper)
        return z::t
    end
end