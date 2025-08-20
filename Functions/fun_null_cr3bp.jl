function fun_null_cr3bp(x0, t0, mu)
    """
    This function is used to investigate specific linear dependencies using the state transition matrix (STM) 
    for a given initial state x0 and time t0. It returns a non-zero vector x that satisfies the equation DF * x = 0, 
    as well as the rank of the matrix DF.
    """
    # Embed the STM (6x6 identity matrix) into the initial state
    Phi0 = Matrix{Float64}(I, 6, 6)
    X0 = vcat(x0, reshape(Phi0, 36))  # 6 + 36 = 42-dimensional

    # Time span
    tspan = (0.0, t0)

    # Solve the ODE
    prob = ODEProblem(fun_stm_cr3bp!, X0, tspan, mu)
    sol = solve(prob, Vern9(), abstol=1e-14, reltol=1e-14)

    # Retrieve the state at the final time
    final_state = sol.u[end]
    X = final_state[1:6]
    Phi = reshape(final_state[7:end], 6, 6)
    dx = zeros(Float64, 6)
    fun_cr3bp!(dx, X, mu, 0.0)
    f_x = dx

    # Construct the matrix DF using selected elements from Phi and f_x
    DF = [
        Phi[2,1] Phi[2,3] Phi[2,5] f_x[2];
        Phi[4,1] Phi[4,3] Phi[4,5] f_x[4];
        Phi[6,1] Phi[6,3] Phi[6,5] f_x[6];
    ]

    # Compute the rank of DF
    R = rank(DF)
    
    # Find a non-trivial vector x that satisfies DF * x = 0
    N = vec(nullspace(DF))
    return N, R
end
