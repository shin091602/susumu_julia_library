function fun_MulShoot_const_CR3BP_cpy(X, r0, rf, n, ToF, mu)
    """
    This is a function used in the Circular Restricted Three-Body Problem (CR3BP) that imposes
    a constraint ensuring that the total flight time from the initial point `r0` to the final point `rf`,
    divided into multiple trajectory segments, matches a predefined Time of Flight (ToF).
    """
    # X       : design variable vector
    # r0      : position of the initial point
    # rf      : position of the final point
    # n       : the number of the patch points
    # ToF     : Time of Flight from the initial to final points
    # mu      : mass ratio of the primaries

    # Equality constraint: total time must equal Time of Flight
    ceq = [sum(X[end - n + 2:end]) - ToF]
    
    # No inequality constraint
    c = Float64[]
    
    return c, ceq
end
