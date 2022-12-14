# Small Prototype for a new tree
using JuMP
import Base.show

# thought: if I have a NodeID, I could have an external map from NodeID to extra info instead of storing it inside the node...?

#  begin with an expression :(max(x,y))
abstract type Node end 
abstract type OpNode <: Node end

mutable struct VarNode <: Node 
    symbol::Symbol
    jump #::JuMP.VariableRef
    val::Float64
end
show(io::IO, a::VarNode) = print(io, "VarNode $(a.symbol)")
mutable struct ConstNode <: Node 
    val::Float64
end
show(io::IO, a::ConstNode) = print(io, "ConstNode $(a.val)")

mutable struct AffOpNode <: OpNode
    # like +, -, *
    symbol::Symbol # operator
    jump #::JuMP.VariableRef
    val::Float64
    children::Array{Node} # ASSUME only 1 or 2 children 
end
show(io::IO, a::OpNode) = print(io, "OpNode $(a.symbol) \n $(a.children)")

# I thought maybe I could do it with an affine node and an OpNode...but I think I need more granular types liek variable and constant
# mutable struct AffNode <: OpNode # afffine expressions 
#     expr::Expr
#     jump::JuMP.AffExpr
#     children::Array{Node}
#     val::Float64
# end

mutable struct PWLOpNode <: OpNode
    symbol::Symbol # operator like max or min or relu
    jump #::JuMP.VariableRef # output variable
    jump_ind #::Array{JuMP.VariableRef} # indicator binary vars 
    val::Float64
    children::Array{Node}
end

"""
A function to replace eval. Hopefully faster.
    If it's not fast enough, I might gain some speed by changing OpNode to more specific node types? And having separate apply_fun functions for them. so I don't have to evaluate this if statement each time. 
"""
function apply_fun(op, args)
    if op == :+
        return +(args...)
    elseif op == :-
        return -(args...)
    elseif op == :*
        return *(args...)
    elseif op == :max
        return max(args...)
    elseif op == :min
        return min(args...)
    elseif op == :relu
        return relu(args...)
    end
end

"""
A function to turn an Expr into a Computational Graph with Nodes.
"""
function populate_graph(e::Expr)
    @assert e.head == :call 
    # args 
    wrapped_args = [populate_graph(a) for a in e.args[2:end]]
    # return parent node
    if e.args[1] ∈ [:+, :-, :*, :/]
        return AffOpNode(e.args[1], "", 0.0, wrapped_args)
    elseif e.args[1] ∈ [:max, :min, :relu]
        return PWLOpNode(e.args[1], "", [], 0.0, wrapped_args)
    else
        warn("Whoa bro wassup here...")
    end 
end
function populate_graph(s::Symbol)
    return VarNode(s, "", 0.0)
end
function populate_graph(r::T where T <: Real)
    return ConstNode(r)
end

"""
A function to evaluate a sample. 
"""
function eval_sample(n::OpNode, sample; mode="custom")
    # evaluate arguments 
    args = [eval_sample(a, sample) for a in n.children]
    if mode == "custom"
        return apply_fun(n.symbol, args)
    elseif mode == "eval"
        return eval(Expr(:call, n.symbol, args...))
    elseif mode == "symengine"
        return Basic(Expr(:call, n.symbol, args...))
    elseif mode == "getproperty"
        return getproperty(Base, n.symbol)(args...)
    end
end 
eval_sample(n::VarNode, sample) = sample[n.symbol]
eval_sample(n::ConstNode, sample) = n.val


# testing ~~~~~~~~~~~~~~~ 

# begin with Expr and population this type
e = :(max(min(x,y),2.45))
g = populate_graph(e)

# to populate this first? And then parse it into a JuMP structure? Or to start with an expr and wrap it into this only as it's parsed into jump? I am leaning towards the second 

# test eval 
x_s, y_s = rand(), rand()
## first with normal functions 
t1 = @elapsed max(min(x_s, y_s), 2.45) 
## then evaluating the graph 
t2 = @elapsed eval_sample(g, Dict(:x => x_s, :y => y_s); mode="custom")
## Compared to running eval in graph
t2p5 = @elapsed eval_sample(g, Dict(:x => x_s, :y => y_s); mode="eval")
## compare to running eval natively 
t3 = @elapsed eval(:(max(eval(:(min(x_s,y_s))),2.45)))
#  compare to running SymEngine 
t4 = @elapsed eval_sample(g, Dict(:x => x_s, :y => y_s); mode="symengine")
# compare to getproperty 
t5 =@elapsed eval_sample(g, Dict(:x => x_s, :y => y_s); mode="getproperty")

println("t1 native = $t1, \n t2 custom = $t2, \n t2p5 eval in graph = $t2p5, \n t3 eval = $t3, \n t4 symengine = $t4, \n t5 getproperty = $t5")
#### It's definitely faster than eval (about 6x) but slower than native code (also about 6x)
# one run: 
# t1 native = 1.1577e-5, t2 custom = 3.3504e-5, t2p5 eval in graph = 0.000177618, t3 eval = 0.000284704, t4 symengine = 9.8378e-5

# Another run, getproperty is super fast!! :D 
"""
t1 native = 1.6258e-5, 
 t2 custom = 8.0403e-5, 
 t2p5 eval in graph = 0.000195858, 
 t3 eval = 0.000318521, 
 t4 symengine = 0.000109191, 
 t5 getproperty = 2.9843e-5
"""

# test encoding into jump
## plan: write 

# test seeding solution 

