function fun_compute_manifolds(L, mu, tf)
    sigma = (1 - mu) / abs(L[1] + mu)^3 + mu / abs(L[1] - 1 + mu)^3
    A = [0 0 0 1 0 0;
         0 0 0 0 1 0;
         0 0 0 0 0 1;
         2 * sigma + 1 0 0 0 2 0;
         0 1 - sigma 0 -2 0 0;
         0 0 -sigma 0 0 0]

    # Obtain eigenvalues and eigenvectors 
    # (Note the order: In MATLAB, `eig(A)` returns `[V, D] = eig(A)`. But in Julia, `D, V = eigen(A)`, where `D` is output as an array.)
    D, V = eigen(A)

    # 
    Vs = zeros(6) 
    Vu = zeros(6) 

    for i in 1:6
        # eigenvalues
        λ = D[i]
        # Extract only real eigenvalues
        if imag(λ) == 0 && real(λ) != 0
            if real(λ) < 0
                # Add eigenvectors corresponding to negative eigenvalues (Convert ComplexF64 → Float64)
                Vs = real(V[:, i])  
            elseif real(λ) > 0
                # Add eigenvectors corresponding to positive eigenvalues (Convert ComplexF64 → Float64)
                Vu = real(V[:, i])  
            end
        end
    end

    # add small perturbation
    x0_sp = vcat(L, zeros(3)) + 1e-10 * Vs / norm(Vs)
    x0_sm = vcat(L, zeros(3)) - 1e-10 * Vs / norm(Vs)
    x0_up = vcat(L, zeros(3)) + 1e-10 * Vu / norm(Vu)
    x0_um = vcat(L, zeros(3)) - 1e-10 * Vu / norm(Vu)

    # Specification integration time 
    tspan_s = (tf, 0)
    tspan_u = (0, tf)

    # Definition of ODEProblem
    prob_sp = ODEProblem(fun_cr3bp!, x0_sp, tspan_s, mu)
    prob_sm = ODEProblem(fun_cr3bp!, x0_sm, tspan_s, mu)
    prob_up = ODEProblem(fun_cr3bp!, x0_up, tspan_u, mu)
    prob_um = ODEProblem(fun_cr3bp!, x0_um, tspan_u, mu)

    # Execution of the ODE Solver
    sol_sp = solve(prob_sp, Vern7(), abstol=1e-14, reltol=1e-14)
    sol_sm = solve(prob_sm, Vern7(), abstol=1e-14, reltol=1e-14)
    sol_up = solve(prob_up, Vern7(), abstol=1e-14, reltol=1e-14)
    sol_um = solve(prob_um, Vern7(), abstol=1e-14, reltol=1e-14)

    return sol_sp, sol_sm, sol_up, sol_um
end