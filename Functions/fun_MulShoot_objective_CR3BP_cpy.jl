function fun_MulShoot_objective_CR3BP_cpy(X, r0, rf, n, mu)
    # X       : design variable vector
    # r0      : position of the initial point
    # rf      : position of the final point
    # n       : the number of the patch points
    # mu      : mass ratio of the primaries

    r0 = vec(r0)
    rf = vec(rf)

    # Reconstruct the full trajectory state vector
    X0 = vcat(r0, X[1:(6n - 9)], rf, X[(6n - 8):end])
    Xf = zeros(6 * (n - 1))

    parameter = (mu)

    # Integrate each trajectory segment
    for i in 1:(n - 1)
        x_init = value.(X0[(6i - 5):(6i)])  # Convert JuMP variables to numerical values
        dt = value(X0[6n + i])             # Convert JuMP variable to numerical value
        prob = ODEProblem(fun_cr3bp!, x_init, (0.0, dt), parameter)  # Pass `mu` as parameter
        sol = solve(prob, abstol=3e-14, reltol=1e-14)
        Xf[(6i - 5):(6i)] = sol[end]
    end

    # Compute mismatch at patch points
    F = zeros(6 * (n - 1))
    for i in 1:(n - 1)
        F[(6i - 5):(6i)] = Xf[(6i - 5):(6i)] - value.(X0[(6i + 1):(6 * (i + 1))])  # Convert JuMP variables to numerical values
    end

    # Objective: norm of mismatches
    return norm(F)
end
