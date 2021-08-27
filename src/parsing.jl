#############################################################
##### Section: Affine #####
#############################################################
function convert_affine_to_jump(s::Symbol, m::Model)
    #println("s: $s")
    if haskey(m, s)
        println("retrieving key: $s using symbol")
        return m[s]
    elseif !isnothing(JuMP.variable_by_name(m, string(s)))
        println("retrieving key: $s using string name")
        return JuMP.variable_by_name(m, string(s))
    else
        println("creating new key: $s")
        # fields  are: (has_lb, lb, has_ub, ub, has_fix, fixed_val, has_start, start, is_binary, is_integer)
        info = VariableInfo(false, NaN, false, NaN, false, NaN, false, NaN, false, false)
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