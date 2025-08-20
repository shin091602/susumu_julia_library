function fun_differential_correction_cr3bp(x0, t0, mu)
    # Initial state and state transition matrix
    X0 = vcat(x0, reshape(I(6), 36))
    tspan = (0.0, t0)
    parameter = (mu)
    
    # Solve the ODE
    prob = ODEProblem(fun_stm_cr3bp!, X0, tspan, parameter)
    sol = solve(prob, Vern7(), reltol=3e-14, abstol=1e-14)
    
    # Extract the final state and state transition matrix
    X = sol.u[end][1:6]

    Phi = reshape(sol.u[end][7:end], 6, 6)
    dx = zeros(6)
    fun_cr3bp!(dx, X, mu, 0.0)
    f_x = dx
    
    # Compute the correction matrix DF
    DF = [Phi[2,3] Phi[2,5] f_x[2];
          Phi[4,3] Phi[4,5] f_x[4];
          Phi[6,3] Phi[6,5] f_x[6]]

    # Compute corrected values
    x_ast = [x0[3]; x0[5]; t0] - DF \ [X[2]; X[4]; X[6]]
    
    # Construct corrected initial state
    x_n = [x0[1]; 0.0; x_ast[1]; 0.0; x_ast[2]; 0.0]
    t_n = x_ast[3]
    
    # Compute Jacobi constant
    C = fun_Jacobi_const(x_n, mu)
    
    return x_n, t_n, C
end