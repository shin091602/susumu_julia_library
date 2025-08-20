function fun_manifold_cr3bp(mu, x0, t0, N, xpert)
    # mu      : mass ratio of the primaries
    # x0      : initial state of a periodic orbit
    # t0      : the period of a periodic orbit
    # N       : the number of points on a periodic orbit
    # xpert   : weight for a position of the eigenvector

    # Define time span
    tspan = range(0, stop=t0, length=N)
    
    # Initial condition including state transition matrix (identity matrix)
    X0 = vcat(x0, reshape(I(6), 36))

    parameter = (mu)

    # Solve ODE
    prob = ODEProblem(fun_stm_cr3bp!, X0, (0.0, t0), parameter)
    sol = solve(prob, Vern7(), saveat=tspan, abstol=3e-14, reltol=1e-14)
    
    # Extract final monodromy matrix
    M = reshape(sol[end][7:42], 6, 6)
    
    # Eigenvalue decomposition
    D, V = eigen(M)
    
    # Find stable and unstable eigenvectors
    index_s = argmin(abs.(D))
    index_u = argmax(abs.(D))
    
    if imag(D[index_s]) != 0 || imag(D[index_u]) != 0
        error("Imaginary eigenvalues are dominant")
    end
    
    vector_stable = V[:, index_s]
    vector_stable *= sign(vector_stable[1])
    
    vector_unstable = V[:, index_u]
    vector_unstable *= sign(vector_unstable[1])
    
    # Initialize manifolds
    XS_left  = zeros(6, N)
    XS_right = zeros(6, N)
    XU_left  = zeros(6, N)
    XU_right = zeros(6, N)
    
    # Compute manifolds
    for iteration in 1:N
        # Extract STM and state at the fixed point
        Phi = reshape(sol[iteration][7:42], 6, 6)
        x_star = sol[iteration][1:6]
        
        # Normalize stable and unstable directions
        S = Phi * vector_stable
        S /= norm(S)
        U = Phi * vector_unstable
        U /= norm(U)
        
        # Apply perturbation
        pert = fill(xpert, 6)
        
        XS_left[:, iteration]  = x_star - S .* pert
        XS_right[:, iteration] = x_star + S .* pert
        XU_left[:, iteration]  = x_star - U .* pert
        XU_right[:, iteration] = x_star + U .* pert
    end
    
    return XS_left, XS_right, XU_left, XU_right, sol
end
