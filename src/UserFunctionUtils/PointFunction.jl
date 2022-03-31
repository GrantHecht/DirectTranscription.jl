# All point functions should be of the form:
#           f!(out, x, u, p, t)
# which looks exactly the same as path functions
# in DirectTranscription.jl. The key difference 
# is PointFunctions can correspond to an arbitrary 
# number of phases with their own unique states,
# controls, and static parameters. Point functions 
# can also involve different instances in time 
# (i.e., the same point function can be evaluated 
# with the state at the initial and final time of 
# a phase).

# When a point function is defined, a set of points
# at which the point function is applied must also 
# be specified. For generality, for a point function
# evaluated at N different points, these points 
# (P_1, P_2, ..., P_N) must be defined. These 
# points can be defined at the initial or final 
# time of any phase in the Trajectory.

# The arguments of the point function x, u, p, and 
# t are defined as:
#
#   x = vcat(x_1, x_2, ..., x_N) where x_j is the 
#   state vector at point P_j
#
#   u = vcat(u_1, u_2, ..., u_N) where u_j is the 
#   control vector at point P_j
#
#   p = vcat(p_1, p_2, ..., p_N) where p_j is the 
#   static parameter vector from the phase 
#   corresponding to P_j
#
#   t = [t_1, t_2, ..., t_N] where t_j is the time
#   corresponding to P_j

abstract type PointFunction end