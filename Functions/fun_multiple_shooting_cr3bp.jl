function fun_multiple_shooting_cr3bp(X0, n, mu)

    Xf = zeros(6*(n-1))
    Phi = Array{Float64, 3}(undef, 6, 6, n-1)

    for i in 1:(n - 1)
        # Embed the STM (6×6 identity matrix) into the initial condition
        Y0 = vcat(X0[(6*i-5):(6*i)], reshape(Matrix(I, 6, 6), :))
        tspan = (0.0, X0[6*n+i])
        # Solve the ODE
        prob = ODEProblem(fun_stm_cr3bp!, Y0, tspan, mu)
        sol = solve(prob, Vern7(), abstol=3e-14, reltol=1e-14)

        Xf[(6*i-5):(6*i)] .= sol.u[end][1:6]
        Phi[:, :, i] = reshape(sol.u[end][7:end], 6, 6)
    end

    # Compute F
    F = zeros(6*(n-1))
    for i in 1:(n - 1)
        F[(6*i-5):(6*i)] .= Xf[(6*i-5):(6*i)] - X0[(6*i+1):(6*(i+1))]
    end

    # Compute DF (the Jacobian matrix)
    DF = zeros(6 * (n - 1), 7 * n - 1)
    I6 = Matrix(I, 6, 6)

    for i in 1:(n-1)
        DF[(6*i-5):(6*i), (6*i-5):(6*i)] .= Phi[:, :, i]
        DF[(6*i-5):(6*i), (6*i+1):(6*(i+1))] .= -I6
        dx = zeros(6)
        fun_cr3bp!(dx, Xf[(6*i-5):(6*i)], mu, 0)
        DF[(6*i-5):(6*i), 6*n+i] .= dx
    end

    # Compute X_ast using Newton-Raphson method
    X_ast = [X0[4:(6*(n-1))]; X0[(6*n-2):end]] - [DF[:,4:6*(n-1)] DF[:,(6*n-2):end]] \ F

    # Reconstruct X_n
    X_n = reshape(vcat(X0[1:3], X_ast[1:(6*n-9)], X0[(6*n-5):(6*n-3)], X_ast[(6*n-8):end]), :, 1)
    
    return X_n
end
