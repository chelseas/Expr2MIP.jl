using SymEngine
using IntervalArithmetic
using LazySets
import JuMP.MOI.OPTIMAL, JuMP.MOI.INFEASIBLE

function Basic2Expr(e)
    return Meta.parse.(string.(expand.(e)))
end

# Apply a symbolic function elementwise to a list of Basics
function elementwise_apply(symbol::Symbol, items::Array)
    # Convert to a list of expressions
    expressions = Basic2Expr(items)
    # Wrap each expression in the function 
    wrapped_expressions = [:($symbol($expr)) for expr in expressions]
    # Convert back to a list of basics
    return Basic.(wrapped_expressions)
end

function apply_multivariate(symbol::Symbol, args::Array)
    # apply a function to multiple args where the args are Basic
    # convert the expression
    expressions = Basic2Expr(args)
    wrapped_expression = Expr(:call, symbol, expressions...)
    return Basic(wrapped_expression)
end

"""
    find_bounds(model, objective)

For a model and objective compute upper and lower bounds on the objective value.
This will typically consist of relaxing a MIP program to an LP and solving that LP.
"""
function find_bounds(model, objective::VariableRef, lp_relaxation=true; bound_type="interval")
    if has_lower_bound(objective) && has_upper_bound(objective)
        @debug "looking up bounds for $objective"
        return lower_bound(objective), upper_bound(objective)
    else
        if bound_type == "opt"
            return find_bounds_through_opt(model, objective, lp_relaxation=lp_relaxation)
        elseif bound_type == "interval"
            return find_bounds_through_int_arithmetic(model, objective)
        end
    end
end

function find_bounds(model, objective::T where T <: Real, lp_relaxation=true; bound_type="interval")
    """
    When you are passed a real value instead of an expression, just return the real value. 
    """
    @debug("For constant objective $objective, returning $objective as up and low bound")
    return objective, objective
end

function find_bounds(model, objective, lp_relaxation=true; bound_type="interval")
    @debug "finding bounds for objective: $objective"
    if bound_type == "opt"
        return find_bounds_through_opt(model, objective, lp_relaxation=lp_relaxation)
    elseif bound_type == "interval"
        return find_bounds_through_int_arithmetic(model, objective)
    end
end

function find_bounds_through_int_arithmetic(model, objective)
    # use interval arithmetic to get bounds on objective
    # turn all variables in expression into intervals (using bounds)
    expr = objective
    int_expr = expr.constant # interval 
    #int_expr_lazy = LazySets.Interval(expr.constant, expr.constant)
    for (variable, coeff) in expr.terms
        var_int = interval(lower_bound(variable), upper_bound(variable))
        #println("var_int = $var_int")
        #var_int_lazy = LazySets.Interval(lower_bound(variable), upper_bound(variable))
        int_expr += coeff*var_int
        #println("int_expr = $int_expr")
        #int_expr_lazy = coeff*LazySets.MinkowskiSum(int_expr_lazy, var_int_lazy)
    end
    #t1 = @elapsed println("int_expr = $int_expr")
    #t2 = @elapsed println("lazy_int_expr = $int_expr_lazy")
    #println("lazy took $t2 to evaluate vs. $t1 for concrete")
    @debug "Found bounds [$(int_expr.lo),$(int_expr.hi)] through interval arithmetic for object $objective"
    return int_expr.lo, int_expr.hi
end

function find_bounds_through_opt(model, objective; lp_relaxation=true)
    undo_relax = nothing
    if lp_relaxation
        @debug("relaxing LP...") 
        undo_relax = relax_integrality(model)
    end

    @objective(model, Min, objective)
    optimize!(model)
    @assert termination_status(model) == OPTIMAL string("Termination Status for computing lower bound was ", termination_status(model))
    lower = objective_value(model)

    @objective(model, Max, objective)
    optimize!(model)
    @assert termination_status(model) == OPTIMAL string("Termination Status for computing upper bound was ", termination_status(model))
    upper = objective_value(model)

    # If we solved the LP relaxation restore the model's integrality constraints
    if undo_relax !== nothing
        undo_relax()
    end

    @debug("found bounds [$(lower), $(upper)] through opt for $objective")

    return lower, upper
end

unit_step(ẑ) = (sign(ẑ)+1)/2