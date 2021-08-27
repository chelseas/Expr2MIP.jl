# this  allows us to dispatch on value of an expression, e.g. we can dispatch on :min
# instead of just being able to dispatch on Symbol or Expr
struct Sym_f{Val} # change to Sym_f{x}
end
Sym_f(x) = Sym_f{x}()
# e.g. Sym_f(:abs)