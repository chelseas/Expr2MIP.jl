using JuMP
import JuMP.MOI.OPTIMAL, JuMP.MOI.INFEASIBLE
using OVERT

include("types.jl")
include("encodings.jl")
include("utilities.jl")

########
# Data Structure to Hold Parameters 
########
mutable struct EncodingParameters
    bound_type::String 
    rel_error_tol::Float64 # for OVERT calls
end

EncodingParameters() = EncodingParameters("interval", 1e-2)

############################################################
#### High Level Functions #####
############################################################

function add_constraint!(model, c::T where T<:Real, var::Symbol; params=EncodingParameters(), expr_map=Dict())
    """
    Adds constraints of the form c == var
    """
    @debug "adding constraint $c::Real"
    c_jumpified = breakdown_and_encode!(model, c, params=params, expr_map=expr_map)

    v = get_or_make_jump_var(model, var, c_jumpified, params.bound_type)

    con_ref = @constraint(model, v == c_jumpified) # actually add the constraint, at last!
    return con_ref, v
end
# the high-level function
function add_constraint!(model, c::Expr, var::Symbol; params=EncodingParameters(), expr_map=Dict())
    """
    Adds constraints of the form c == var
    Every symbol in the expression c must exist, but the variable var may not exist yet.
    """
    @debug "adding constraint $var == $c::Expr"
    println("encoding parameters are: $params")
    # adds constraints to the jump model of the form: var == c
    # where c may be an arbitrary expression like max(10*u, step(z)*u + 6) - 12
    c = convert_step_times_var(c)
    c_jumpified = breakdown_and_encode!(model, c, params=params, expr_map=expr_map)
    # c_jumpified may come back as: v_46 or v_12 - 25 + 14
    
    v = get_or_make_jump_var(model, var, c_jumpified, params.bound_type)

    #println("v = $v")
    #println("typeof(v)= $(typeof(v))")
    con_ref = @constraint(model, v == c_jumpified) # actually add the constraint, at last!
    return con_ref, v
end

function add_constraint!(model, c::Symbol, var::Symbol; params=EncodingParameters(), expr_map=Dict())
    """
    Add constraints of the form c::Symbol == var.
    c must exist in the model already, var may not exist yet.
    """
    @debug "adding constraint $var == $c::Symbol"
    c_jumpified = breakdown_and_encode!(model, c, params=params, expr_map=expr_map)
    v = get_or_make_jump_var(model, var, c_jumpified, params.bound_type)
    con_ref = @constraint(model, v == c_jumpified) # actually add the constraint, at last!
    return con_ref, v
end

# helper 
# TODO: this function is much the same as convert_affine_to_jump for var. could they be combined in some way?
function get_or_make_jump_var(model, var, c_jumpified, bound_type)
    if haskey(model, var)
        @debug "retrieving symbol key $var"
        v = model[var]
    elseif !isnothing(JuMP.variable_by_name(model, string(var)))
        @debug "retrieving string key $var"
        v = JuMP.variable_by_name(model, string(var))
        # TODO: one thing to do that would be interesting would be to run the LP relaxation and see if the bounds are tighter than the OVERT ones...
    else # variable does not exist yet
        # Find upper and lower bounds then create the variable
        @debug "creating new variable $var"
        lower, upper = find_bounds(model, c_jumpified, bound_type=bound_type)
        v = @variable(model, base_name=string(var), lower_bound=lower, upper_bound=upper) # output variable
    end
    return v
end

# max(x + y + 7z, y + z)
function breakdown_and_encode!(model, expr::Expr; params=EncodingParameters(), expr_map=Dict())
    """
    Main Parser.
    This is the main function that gets called inside the encode! functions in order to parse arguments. 
    """
    @debug "breakdown and encode $expr::Expr"
    #containers = (constraints, output_variables, bounds, model)
    #println("expr: $expr")

    # If the expr is already present in the model, just return the jump reference
    if haskey(expr_map, expr)
        return expr_map[expr]
    end

    f = expr.args[1]
    args = expr.args[2:end]
    # base cases
    if is_affine(expr)
        @debug "is affine"
        # encode
        jumped_expr = convert_affine_to_jump(expr, model)
        #println("jumped_expr=$(jumped_expr)")
        # add mapping from expression expr to jump reference 
        expr_map[expr] = jumped_expr
        return jumped_expr
    elseif f ∈ [:abs, :max, :min, :relu, :unit_step, :unit_step_times_var, :+, :-, :*, :/]
        @debug "encoding PWL or affine function $f"
        out = encode!(model, Sym_f(f), args; params=params, expr_map=expr_map) # returns jump compatible type
        # add mapping from expression expr to jump reference 
        expr_map[expr] = out
        return out
    else # assuming that this should be passed to OVERT
        println("calling OVERT for expr $expr")
        # smooth nonlinearity
        out = call_overt!(model, f, args; params=params, expr_map=expr_map)
        # add mapping from expression expr to jump reference 
        expr_map[expr] = out
        return out
    end
end

function breakdown_and_encode!(model, s::Union{Symbol, T where T <: Real}; params=EncodingParameters(), expr_map=Dict())
    @debug "breakdown and encode symbol or number $s"
    return convert_affine_to_jump(s, model)
end

############################################################
###### recursive cases ######
############################################################
##################################################################################
# Encode into the jump model ####################################################
##################################################################################
# GenericAffineExpr -> JumpVariableRef dict for memoization later?
function encode!(model, wrapped_f::Sym_f{:abs}, args::Array; params=EncodingParameters(), expr_map=Dict()) # -> JuMP variable ref
    """
    This function handles recursive parsing and encoding of absolute value. 
    """
    @debug "encoding abs of $(args[1])"
    @assert length(args) == 1 # scalar function applied to scalar input
    encoded_args = [breakdown_and_encode!(model, a; params=params, expr_map=expr_map) for a in args]
    input_arg = encoded_args[1]

    lower, upper = find_bounds(model, input_arg, bound_type=params.bound_type)
    output_var = encode_abs!(model, input_arg, lower, upper)
    return output_var
end
function encode!(model, wrapped_f::Sym_f{:max}, args::Array; params=EncodingParameters(), expr_map=Dict(), relu_rewrite=false)
    if !relu_rewrite 
        return encode_max_direct!(model, args::Array; params=params, expr_map=expr_map)
    else # relu_rewrite == true
        return encode_max_relu_rewrite!(model, args::Array; params=params, expr_map=expr_map)
    end
end
function encode_max_direct!(model, args::Array; params=EncodingParameters(), expr_map=Dict())
    @debug "encoding max of $(args) EXACTLY"
    encoded_args = [breakdown_and_encode!(model, a, params=params, expr_map=expr_map) for a in args]
    n_args = length(encoded_args)

    bounds = [find_bounds(model, encoded_arg, bound_type=params.bound_type) for encoded_arg in encoded_args]
    lower_bounds = [bound_tuple[1] for bound_tuple in bounds]
    upper_bounds = [bound_tuple[2] for bound_tuple in bounds]
    y = encode_max_real!(model, encoded_args, lower_bounds, upper_bounds)
    return y
end

function encode_max_relu_rewrite!(model, args::Array; params=EncodingParameters(), expr_map=Dict())
    # use relu-re-write 
    # max(x,y) = 0.5*(x + y + relu(x-y) + relu(y-x))
    # TODO: Why does using the rewrite + triangle relaxation add so many more constraints?

    # trying something to see if it adds speed by reducing the number of constraints. 
    encoded_args = [breakdown_and_encode!(model, a, params=params, expr_map=expr_map) for a in args]
    bounds = [find_bounds(model, encoded_arg, bound_type=params.bound_type) for encoded_arg in encoded_args]
    lower_bounds = [bound_tuple[1] for bound_tuple in bounds]
    upper_bounds = [bound_tuple[2] for bound_tuple in bounds]
    need_max, indices_to_keep, l_max = check_if_need_max(lower_bounds, upper_bounds)

    if !need_max
        @debug "There is no need to actually take max." encoded_args[indices_to_keep]
        @assert length(indices_to_keep) == 1
        return encoded_args[indices_to_keep...]
    else 
        @assert length(args) == 2 
        @debug "Encoding max of $(args) using relu-rewrite"
        x = args[1]
        y = args[2]
        # TODO: could probably re-use some terms and/or some bounds because (x-y) == -(y-x).
        # Maybe could add a common term? to connect relu(x-y) and relu(y-x) ?
        new_expr = :(0.5*($x + $y + relu($x-$y) + relu($y-$x)))
        return breakdown_and_encode!(model, new_expr, params=params, expr_map=expr_map)
    end
end
function encode!(model, wrapped_f::Sym_f{:relu}, args::Array; params=EncodingParameters(), expr_map=Dict())
    @debug "Encoding relu of $(args)"
    encoded_args = [breakdown_and_encode!(model, a, params=params, expr_map=expr_map) for a in args]
    @assert length(encoded_args) == 1 # relu is a 1-arg function 
    lower, upper = find_bounds(model, encoded_args[1], bound_type=params.bound_type)

    y = encode_relu!(model, encoded_args[1], lower, upper)
    return y
end
function encode!(model, wrapped_f::Sym_f{:min}, args::Array; params=EncodingParameters(), expr_map=Dict())
    @debug "encoding min of $(args)"
    return -encode!(model, Sym_f(:max), [:(-$a) for a in args], params=params, expr_map=expr_map)
end
function encode!(model, wrapped_f::Sym_f{:unit_step}, args::Array; params=EncodingParameters(), expr_map=Dict())
    @debug "encoding unit step of $(args)"
    encoded_args = [breakdown_and_encode!(model, a, params=params, expr_map=expr_map) for a in args]
    @assert length(encoded_args) == 1

    lower, upper = find_bounds(model, encoded_args[1], bound_type=params.bound_type)
    δ = encode_unit_step!(model, encoded_args[1], lower, upper)
    return δ 
end
function encode!(model, wrapped_f::Sym_f{:unit_step_times_var}, args::Array; params=EncodingParameters(), expr_map=Dict())
    @debug "encode unit step times var of $(args)"
    ẑ, x = args
    # z = unit_step(ẑ)*x
    ẑₑ = breakdown_and_encode!(model, ẑ, params=params, expr_map=expr_map)
    xₑ = breakdown_and_encode!(model, x, params=params, expr_map=expr_map)
    δ = @variable(model, binary=true, base_name="δ_ustv")

    # TODO: @castrong add computation of lower and upper bounds of both input to the unit step and the real valued variable 
    lower_ẑₑ, upper_ẑₑ = find_bounds(model, ẑₑ, bound_type=params.bound_type)
    lower_xₑ, upper_xₑ = find_bounds(model, xₑ, bound_type=params.bound_type)
    z = encode_unit_step_times_var!(model, ẑₑ, xₑ, δ, lower_ẑₑ, upper_ẑₑ, lower_xₑ, upper_xₑ)
    return z
end
# version that logs bounds in a dict?
# function encode!(model, wrapped_f::Sym_f{:unit_step_times_var}, args::Array, ranges::Dict, δ)
#     ẑ, x = args
#     # z = unit_step(ẑ)*x
#     ẑₑ = breakdown_and_encode!(model, ẑ, ranges)
#     xₑ = breakdown_and_encode!(model, x, ranges)
#     l, u = ranges[ẑ]
#     γ, ζ = ranges[x]
#     z = encode_unit_step_times_var!(model, ẑₑ, xₑ, δ, l, u, γ, ζ)
#     return z
# end
function encode!(model, wrapped_f::Sym_f{:+}, args::Array; params=EncodingParameters(), expr_map=Dict())
    encoded_args = [breakdown_and_encode!(model, a, params=params, expr_map=expr_map) for a in args]
    return +(encoded_args...)
end
function encode!(model, wrapped_f::Sym_f{:-}, args::Array; params=EncodingParameters(), expr_map=Dict())
    encoded_args = [breakdown_and_encode!(model, a, params=params, expr_map=expr_map) for a in args]
    return -(encoded_args...)
end
function encode!(model, wrapped_f::Sym_f{:*}, args::Array; params=EncodingParameters(), expr_map=Dict())
    println("Dealing with multiplication of args: $(args)")
    isnn = .!is_number.(args)
    if sum(isnn) > 1 # multiplication of two real valued variables is present
        # todo: check we don't have e.g. relu*relu or abs*abs...because I still have to add support to overt for relu*relu
        # also check for unit_step*unit_step and steptimesvar*steptimesvar
        var = call_overt!(model, :*, args[isnn], params=params, expr_map=expr_map) # encodes this arg(s) that contain variables. returns jump var
        return *(var, args[.!isnn]...) # multiply variable and coefficient together
    else # "outer affine" e.g. const*relu
        encoded_args = [breakdown_and_encode!(model, a, params=params, expr_map=expr_map) for a in args[isnn]]
        return *(encoded_args..., args[.!isnn]...)
    end
end
function encode!(model, wrapped_f::Sym_f{:/}, args::Array; params=EncodingParameters(), expr_map=Dict())
    @assert length(args) == 2
    if !is_number(args[2]) # call overt to handle c/x or c/(x+y) type stuff
        var = call_overt!(model, :/, args; params=params, expr_map=expr_map)
        return var
    else  # looks something like: relu(x) / c, perhaps. "outer affine"
        # is_number(args[2]) # second arg is number
        arg1 = breakdown_and_encode!(model, args[1], params=params, expr_map=expr_map)
        return arg1/args[2] # use operator overloadding to construct GenericAffExpr
    end   
end

#############################################################
##### Section: Affine #####
#############################################################
special_symbols = Set{Symbol}([:pi])

function convert_affine_to_jump(s::Symbol, m::Model)
    @debug("converting affine symbol: $s to JuMP")
    if s in special_symbols
        @debug("processing special symbol $s")
        return eval(s)
    elseif haskey(m, s)
        @debug("retrieving key: $s using symbol")
        return m[s]
    elseif !isnothing(JuMP.variable_by_name(m, string(s)))
        @debug("retrieving key: $s using string name")
        return JuMP.variable_by_name(m, string(s))
    else
        @warn "don't think we should be using this case...investigate!"
        @assert 1 == 0
        @debug("creating new key: $s")
        # fields  are: (has_lb, lb, has_ub, ub, has_fix, fixed_val, has_start, start, is_binary, is_integer)
        info = VariableInfo(false, NaN, false, NaN, false, NaN, false, NaN, false, false)

        # TODO: @castrong add bounds in here 
        JuMP.add_variable(m, JuMP.build_variable(error, info), string(s))
        #return @variable(m, $s) # then m[s] should return the corresponding JuMP variable
    end
end
function convert_affine_to_jump(c::T where T <: Real, m::Model)
    #println("c: $c")
    return c
end 
# todo: worry about log(3) coefficients -> JK! JuMP can handle these! :)
function convert_affine_to_jump(expr::Expr, m::Model)
    @debug "converting affine expr: $expr to JuMP"
    #println("Expr: $expr")
    f = expr.args[1]
    args = expr.args[2:end]
    jump_args = [convert_affine_to_jump(a, m) for a in args]
    #println("f= $f")
    #println("args= $args")
    return eval(f)(jump_args...) # return GenericAffExpr or other JuMP compatible type
end
###############################################################

# codee swiped from OVERT Overapprox repo ---------------------------
function is_affine(n::T where T <: Real)
    return true
end

function is_affine(s::Symbol)
    return true
end

function is_number(expr)
    try
        e = eval(expr)
        return e isa T where T <: Real
    catch
        return false
    end
end

function is_affine(expr::Expr)
    """
    given an expression expr, this function determines if the expression
    is an affine function.

    Example: is_affine(:(x+1))       = true
             is_affine(:(-x+(2y-z))) = true
             is_affine(:(log(x)))    = false
             is_affine(:(x + x*z))   = false
             is_affine(:(x/6))       = true
             is_affine(:(5*x))       = true
             is_affine(:(log(2)*x))  = true
             is_affine(:(-x))        = true
     """
    # it is number
    if is_number(expr)
        @debug "$expr is number "
        return true
    else
        func = expr.args[1]
        args = expr.args[2:end]
        if !(func ∈ [:+, :-, :*, :/]) # only these operations are allowed
            @debug "$func not in allowed set"
            return false
        else  # func ∈ [:+, :-, :*, :/]
            @debug "$func ∈ [:+, :-, :*, :/]"
            if func == :* # one arg can be affine and rest must be a number
                @debug "testing if mul expr $expr is affine"
                options = [is_affine(args[i]) && all(is_number.(args[[1:i-1...,i+1:end...]])) for i in 1:length(args)]
                return any(options)
            elseif func == :/ # only two args allowed. second arg has to be a number
                @debug "testing if div expr $expr is affine"
                return is_affine(args[1]) && is_number(args[2])
            else # func is + or -
                @debug "testing if +/- expr $expr is affine"
                return all(is_affine.(args))
            end
        end
    end
end

# -------------------------------------------------------------------------------

############################################################
#####  overt interface ######
############################################################
function define_state_variables!(model, domain)
    for (state,bounds) in domain 
        if !has_key(model, state)
            @variable(model, base_name=string(state), lower_bound=bounds[1], upper_bound=bounds[2])
        else
            @debug("skipping. don't want to redefine $state ∈ $bounds")
        end
    end
    return nothing # purely a side effect function
end

function encode_overapprox!(model::Model, oa::OverApproximation, domain; params=EncodingParameters(), expr_map=Dict())
    # From the overapproximation object we want to encode the approx_eq and the approx_ineq
    # we need all variables to be defined before they are encoded bc we expect ranges
    println("Encoding the overapproximation using $(params.bound_type) bounding.")
    define_state_variables!(model, domain)
    # approx_eq
    for c in oa.approx_eq
        # constraints are of the form: v_1 == 5*v_6 - 12*v_9 - 7
        LHS = c.args[2] # just a symbol
        @assert typeof(LHS) == Symbol 
        RHS = c.args[3] # an expr
        con_ref, output_ref = add_constraint!(model, RHS, LHS; params=params, expr_map=expr_map)
    end
    # approx_ineq
    for c in oa.approx_ineq
        # these constraints are of the form: v_1 <= v_3
        @assert typeof(c.args[2]) == Symbol 
        @assert typeof(c.args[3]) == Symbol
        LHS = JuMP.variable_by_name(model, string(c.args[2]))
        RHS = JuMP.variable_by_name(model, string(c.args[3]))
        @constraint(model, LHS <= RHS)
    end
    # TODO: use both here and in fuel gage. (don't forget to push to cloud and then pull on jodhpur)
end

function replace_arg(arg::GenericAffExpr)
    global NEW_OVERT_VAR_COUNT += 1
    return Symbol("nov_$(NEW_OVERT_VAR_COUNT)") 
end
replace_arg(arg::VariableRef) = Symbol(name(arg)) 
replace_arg(arg::t where t <: Real) = arg # identity

function replace_arg_in_model(arg::Symbol, lb, ub)
    if has_key(model, string(arg))
        return variable_by_name(model, string(arg))
    else
        var = @variable(model, base_name=string(arg), lower_bound=lb, upper_bound=ub)

    end
end

function call_overt!(model, f, args; params=EncodingParameters(), expr_map=Dict(), N=-1)
    println("Encoding $f applied to $args using OVERT! :) ")
    @debug "Encoding using OVERT"
    # encode args 
    encoded_args = [breakdown_and_encode!(model, a, params=params, expr_map=expr_map) for a in args]
    @debug "Encoded args are " encoded_args 
    # replace args with new variables IF NEEDED
    new_args = [replace_arg(a) for a in encoded_args]
    @debug "New args are " new_args 
    # global NEW_OVERT_VAR_COUNT += length(args) 
    expr = Expr(:call, f, new_args...)
    # find ranges for new args 
    bounds = [[find_bounds(model, encoded_arg, bound_type=params.bound_type)...] for encoded_arg in encoded_args]
    @debug "Bounds for new args are: " bounds 
    # define new args in model IF NEEDED
    for (i, new_arg) in enumerate(new_args)
        if !has_key(model, string(new_arg))   # if it doesn't exist in the model yet (Bc it came from an expression) 
            new_arg_ref = @variable(model, base_name=string(new_arg), lower_bound=bounds[i][1], upper_bound=bounds[i][2])
            @debug "New arg ref is " new_arg_ref
            # connect new arg refs to args in the model that they are equivalent to!!!
            @debug "Adding constraint to make $(encoded_args[i]) == $(new_arg_ref)"
            @constraint(model, encoded_args[i] == new_arg_ref)
        end
    end
    # construct range_dict
    range_dict = Dict(zip(new_args, bounds)) 
    # overapproximate
    @debug "Overapproximating $expr over domain $(range_dict)"
    oa = overapprox(expr, range_dict::Dict{Symbol, Array{T, 1}} where {T <: Real}, N=N, rel_error_tol=params.rel_error_tol)
    @debug "Output variable of overapprox is: $(oa.output) with range $(oa.output_range)"
    # TODO: does this output variable have bounds? Or should I add e.g. those from OVERT?
    # deal with overt
    encode_overapprox!(model, oa, oa.ranges; params=params, expr_map=expr_map)
    @assert has_key(model, string(oa.output))
    return JuMP.variable_by_name(model, string(oa.output)) # which should already be in the model
end

############################################################
##### unit step times var stuff #####
############################################################

function convert_step_times_var(s::Symbol)
    return s
end
function convert_step_times_var(n::T where T <: Real)
    return n
end

function is_a_step_function(s::Symbol)
    return s == :(unit_step)
end
function is_a_step_function(n::T where T <: Real)
    return false
end
function is_a_step_function(e::Expr)
    f = e.args[1]
    return f == :(unit_step)
end

function convert_step_times_var(expr::Expr)
    # any place where the expression unit_step(z)*x appears, replace it with: unit_step_times_var(z, x)
    # recursive: if this is a multiplication, check if any of the terms contain step*var
    @debug "convert step times var"
    f = expr.args[1]
    args = expr.args[2:end]
    if f == :*
        is_step = [is_a_step_function(a) for a in args]
        is_var = [!is_number(a) && is_affine(a) for a in args] # e.g. x or (x + y)
        has_step_inside = [occursin("unit_step", string(a)) for a in args]
        if sum(is_step) > 1 || sum(has_step_inside) > 1
            @error "Multiplication of two unit steps not yet supported"
        elseif sum(is_step) == 1
            step_f_expr = args[is_step][1]
            step_arg = step_f_expr.args[2]
            @assert length(step_f_expr.args[2:end]) == 1 # should only be 1 arg to a step function
            # recurse on args to step function
            step_arg_converted = convert_step_times_var(step_arg)
            vars = args[is_var]
            # recurse on vars
            vars_converted = [convert_step_times_var(v) for v in vars]
            if length(vars) > 1
                var_expr = Expr(:call, :*, vars_converted...)
            else
                var_expr = vars_converted[1]
            end
            unit_step_expr = :(unit_step_times_var($(step_arg_converted), $(var_expr)))
            if length(args) - sum(is_step) - sum(is_var) > 0
                return Expr(:call, :*, args[(.!is_step) .&  (.!is_var)]..., unit_step_expr)
            else
                return unit_step_expr
            end
        elseif sum(is_step) == 0
            converted_args = [convert_step_times_var(a) for a in args]
            return Expr(:call, f, converted_args...)
        end
    else
        converted_args = [convert_step_times_var(a) for a in args]
        return Expr(:call, f, converted_args...)
    end
end

##################################################################
### Helper Functions
##################################################################

function has_key(model::Model, key::Symbol)
    return haskey(model, key) || !isnothing(JuMP.variable_by_name(model, string(key))) 
end
function has_key(model::Model, key::String)
    return !isnothing(JuMP.variable_by_name(model, key)) 
end

