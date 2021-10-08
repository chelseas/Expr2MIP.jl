module Expr2MIP

include("parsing.jl")

export add_constraint!,
       breakdown_and_encode!,
       encode!,
       convert_affine_to_jump,
       is_affine,
       encode_abs!,
       encode_max_real!,
       encode_unit_step!,
       encode_unit_step_times_var!,
       elementwise_apply,
       apply_multivariate,
       Basic2Expr,
       find_bounds
end