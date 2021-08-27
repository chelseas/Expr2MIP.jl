using SymEngine

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