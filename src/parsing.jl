using JuMP
import JuMP.MOI.OPTIMAL, JuMP.MOI.INFEASIBLE

include("types.jl")
include("encodings.jl")
include("utilities.jl")
ENV["JULIA_DEBUG"] = Main

############################################################
#### High Level Functions #####
############################################################

function add_constraint!(model, c::T where T<:Real, var::Symbol)
    """
    Adds constraints of the form c == var
    """
    @debug "adding constraint $c::Real"
    c_jumpified = breakdown_and_encode!(model, c)

    if haskey(model, var)
        @debug "retrieving symbol key $var"
        v = model[var]
    elseif !isnothing(JuMP.variable_by_name(model, string(var)))
        @debug "retrieving string key $var"
        v = JuMP.variable_by_name(model, string(var))
    else # variable does not exist yet
        # Find upper and lower bounds then create the variable
        @debug "creating new variable"
        lower, upper = find_bounds(model, c_jumpified)
        v = @variable(model, base_name=string(var), lower_bound=lower, upper_bound=upper) # output variable
    end

    con_ref = @constraint(model, v == c_jumpified) # actually add the constraint, at last!
    return con_ref, v
end
# the high-level function
function add_constraint!(model, c::Expr, var::Symbol)
    """
    Adds constraints of the form c == var
    """
    @debug "adding constraint $var == $c::Expr"
    # adds constraints to the jump model of the form: var == c
    # where c may be an arbitrary expression like max(10*u, step(z)*u + 6) - 12
    c = convert_step_times_var(c)
    c_jumpified = breakdown_and_encode!(model, c)
    # c_jumpified may come back as: v_46 or v_12 - 25 + 14
    
    if haskey(model, var)
        @debug "retrieving symbol key $var"
        v = model[var]
    elseif !isnothing(JuMP.variable_by_name(model, string(var)))
        @debug "retrieving string key $var"
        v = JuMP.variable_by_name(model, string(var))
        # TODO: one thing to do that would be interesting would be to run the LP relaxation and see if the bounds are tighter than the OVERT ones...
    else # variable does not exist yet
        # Find upper and lower bounds then create the variable
        @debug "creating new variable"
        lower, upper = find_bounds(model, c_jumpified)
        v = @variable(model, base_name=string(var), lower_bound=lower, upper_bound=upper) # output variable
    end

    #println("v = $v")
    #println("typeof(v)= $(typeof(v))")
    con_ref = @constraint(model, v == c_jumpified) # actually add the constraint, at last!
    return con_ref, v
end

# max(x + y + 7z, y + z)
function breakdown_and_encode!(model, expr::Expr)
    @debug "breakdown and encode $expr::Expr"
    #containers = (constraints, output_variables, bounds, model)
    #println("expr: $expr")
    f = expr.args[1]
    args = expr.args[2:end]
    # base cases
    if is_affine(expr)
        @debug "is affine"
        # encode
        jumped_expr = convert_affine_to_jump(expr, model)
        #println("jumped_expr=$(jumped_expr)")
        return jumped_expr
    elseif f ∈ [:abs, :max, :min, :relu, :unit_step, :unit_step_times_var, :+, :-, :*, :/]
        @debug "encoding PWL or affine function $f"
        out = encode!(model, Sym_f(f), args) # returns jump compatible type
        return out
    else # assuming that this should be passed to OVERT
        println("calling OVERT for expr $expr")
        # smooth nonlinearity
        out = call_overt(model, f, args)
        return out
    end
end

function breakdown_and_encode!(model, s::Union{Symbol, T where T <: Real})
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
function encode!(model, wrapped_f::Sym_f{:abs}, args::Array) # -> JuMP variable ref
    """
    This function handles recursive parsing and encoding of absolute value. 
    """
    @debug "encoding abs of $(args[1])"
    @assert length(args) == 1 # scalar function applied to scalar input
    encoded_args = [breakdown_and_encode!(model, a) for a in args]
    input_arg = encoded_args[1]

    lower, upper = find_bounds(model, input_arg)
    output_var = encode_abs!(model, input_arg, lower, upper)
    return output_var
end
function encode!(model, wrapped_f::Sym_f{:max}, args::Array)
    @debug "encoding max of $(args)"
    encoded_args = [breakdown_and_encode!(model, a) for a in args]
    n_args = length(encoded_args)

    bounds = [find_bounds(model, encoded_arg) for encoded_arg in encoded_args]
    lower_bounds = [bound_tuple[1] for bound_tuple in bounds]
    upper_bounds = [bound_tuple[2] for bound_tuple in bounds]
    y = encode_max_real!(model, encoded_args, lower_bounds, upper_bounds)
    return y
end
function encode!(model, wrapped_f::Sym_f{:min}, args::Array)
    @debug "encoding min of $(args)"
    return -encode!(model, Sym_f(:max), [:(-$a) for a in args])
end
function encode!(model, wrapped_f::Sym_f{:unit_step}, args::Array)
    @debug "encoding unit step of $(args)"
    encoded_args = [breakdown_and_encode!(model, a) for a in args]
    @assert length(encoded_args) == 1

    
    lower, upper = find_bounds(model, encoded_args[1])
    δ = encode_unit_step!(model, encoded_args[1], lower, upper)
    return δ 
end
function encode!(model, wrapped_f::Sym_f{:unit_step_times_var}, args::Array)
    @debug "encode unit step times var of $(args)"
    ẑ, x = args
    # z = unit_step(ẑ)*x
    ẑₑ = breakdown_and_encode!(model, ẑ)
    xₑ = breakdown_and_encode!(model, x)
    δ = @variable(model, binary=true, base_name="δ_ustv")

    # TODO: @castrong add computation of lower and upper bounds of both input to the unit step and the real valued variable 
    lower_ẑₑ, upper_ẑₑ = find_bounds(model, ẑₑ)
    lower_xₑ, upper_xₑ = find_bounds(model, xₑ)
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
function encode!(model, wrapped_f::Sym_f{:+}, args::Array)
    encoded_args = [breakdown_and_encode!(model, a) for a in args]
    return +(encoded_args...)
end
function encode!(model, wrapped_f::Sym_f{:-}, args::Array)
    encoded_args = [breakdown_and_encode!(model, a) for a in args]
    return -(encoded_args...)
end
function encode!(model, wrapped_f::Sym_f{:*}, args::Array)
    #println("Dealing with multiplication")
    isnn = .!is_number.(args)
    if sum(isnn) > 1 # multiplication of two real valued variables is present
        # todo: check we don't have e.g. relu*relu or abs*abs...because I still have to add support to overt for relu*relu
        # also check for unit_step*unit_step and steptimesvar*steptimesvar
        var = call_overt(model, :*, args[isnn]) # encodes this arg(s) that contain variables. returns jump var
        return *(var, args[.!isnn]...) # multiply variable and coefficient together
    else # "outer affine" e.g. const*relu
        encoded_args = [breakdown_and_encode!(model, a) for a in args[isnn]]
        return *(encoded_args..., args[.!isnn]...)
    end
end
function encode!(model, wrapped_f::Sym_f{:/}, args::Array)
    @assert length(args) == 2
    if !is_number(args[2]) # call overt to handle c/x or c/(x+y) type stuff
        var = call_overt(model, :/, args)
        return var
    else  # looks something like: relu(x) / c, perhaps. "outer affine"
        # is_number(args[2]) # second arg is number
        arg1 = breakdown_and_encode!(model, args[1])
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
        return true
    else
        func = expr.args[1]
        args = expr.args[2:end]
        if !(func ∈ [:+, :-, :*, :/]) # only these operations are allowed
            return false
        else  # func ∈ [:+, :-, :*, :/]
            if func == :* # one arg can be affine and rest must be a number
                options = [is_affine(args[i]) && all(is_number.(args[[1:i-1...,i+1:end...]])) for i in 1:length(args)]
                return any(options)
            elseif func == :/ # only two args allowed. second arg has to be a number
                return is_affine(args[1]) && is_number(args[2])
            else # func is + or -
                return all(is_affine.(args))
            end
        end
    end
end

# -------------------------------------------------------------------------------

############################################################
#####  overt interface ######
############################################################
function call_overt(model, f, args)
    # todo: turn from sketch -> real code
    expr = Expr(:call, f, args...)
    oa = overapprox_nd(expr, range_dict::Dict{Symbol, Array{T, 1}} where {T <: Real}, N=-1)
    # deal with overt
    for constraint in oa.approx_eq #, oa.approx_ineq
        add_constraint!(model, constraint, :c)
    end
    return oa.output
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

