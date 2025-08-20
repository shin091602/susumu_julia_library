function fun_differential_correction_cr3bp_PAC(x01_ast, t01_ast, x02, t02, scale, delta, mu)
    # Initial state and state transition matrix
    X0 = vcat(x02, reshape(I(6), 36))
    tspan = (0.0, t02)
    parameter = (mu)
    
    # Solve the ODE
    prob = ODEProblem(fun_stm_cr3bp!, X0, tspan, parameter)
    sol = solve(prob, Vern9(), reltol=1e-14, abstol=1e-14)
    
    # Extract the final state and state transition matrix
    X = sol.u[end][1:6]

    Phi = reshape(sol.u[end][7:end], 6, 6)
    dx = zeros(6)
    fun_cr3bp!(dx, X, mu, 0.0)
    f_x = dx
    F = [X[2]; X[4]; X[6]]
    F_PAC = dot(([x02[1]; x02[3]; x02[5]; t02] .- [x01_ast[1]; x01_ast[3]; x01_ast[5]; t01_ast]),delta) - scale
    G = [F; F_PAC]

    # Compute the correction matrix DF
    DF = [  Phi[2, 1]  Phi[2, 3]  Phi[2, 5]  f_x[2];
            Phi[4, 1]  Phi[4, 3]  Phi[4, 5]  f_x[4];
            Phi[6, 1]  Phi[6, 3]  Phi[6, 5]  f_x[6]]

    DG = [DF; delta']
    
    # Compute corrected values
    X_ast = [x02[1]; x02[3]; x02[5]; t02] .- DG\G
    
    # Construct corrected initial state
    x02_ast = [X_ast[1]; 0; X_ast[2]; 0; X_ast[3]; 0]
    t02_ast = X_ast[4]
    C = fun_Jacobi_const(x02_ast, mu)
    
    return x02_ast, t02_ast, C, G 
end