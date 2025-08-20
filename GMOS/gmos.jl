function center_tangent_qpo(zpo, gp)
    # =================================
    # =================================
    #
    # Center Tangent QPO
    # By: Damennick Henry
    # Date: 11/10/21
    # 
    # Description: Compute a linear approximation of a qpo emanating from a po from its center tangent vector
    #
    # Inputs
    #       zpo - Periodic orbit solution
    #       par - Parameters dict
    #
    # Outputs
    #
    # Dependencies:
    #
    # Notes: See  "Fully Numerical Methods for Continuing Families of Quasi-Periodic Invariant Tori in Astrodyanamics" by N. Baresi et al for description of the GMOS algorithm
    #      
    # =================================
    # =================================

    # INITIALIZATION
    # Torus solution dictionary
    Ud0 = Dict();

    # Dictionary extraction
    # State dimension
    d = gp["d"];
    # State and parameter dimension
    dp = gp["dp"];
    # Number of invariant curve points
    N = gp["N"];
    # Number of multiple shooting points
    M = gp["M"];
    # Fourier matrices
    Fr = gp["Fr"]; IFr = gp["IFr"]; DFr = gp["DFr"];

    # Periodic orbit solution extraction
    # Periodic orbit state
    xpo = zpo[1:d];
    # Orbit period
    Tpo = zpo[7];
    # System parameters
    p = zpo[7+1:end];
    # Time vector
    tk = collect(LinRange(0 , Tpo, M + 1)); pop!(tk);
    # Invariant curve angles
    thtj = collect(LinRange(0, 2*pi, N + 1)); pop!(thtj);

    # MONODROMY MATRIX
    # Initial STM
    phi0 = 1.0*Matrix(I, dp, dp);
    # Create ODE problem
    prob = ODEProblem(gp["eom"], [xpo; phi0[:]], (0.0, Tpo), p);
    # Solve ODE
    sol = solve(prob, VCABM(), abstol = 1e-16, reltol = 3e-14);
    #=
    anim = @animate for i = 1:2:length(x)
        plot(x[1:i], y[1:i], legend=false)
    end
    gif(anim, "LPO.gif", fps = 30)
    =#
    #display(plot3d!([u[1] for u in sol.u], [u[2] for u in sol.u], [u[3] for u in sol.u]))
    # Final state
    xT = sol.u[end]
    # Monodromy matrix
    phiT = reshape(xT[d+1:end], (dp, dp));
    phiT = phiT[1:d, 1:d]
    # Eigendata of monodromy matrix
    Me = eigen(phiT);
    # Eigenvalue index corresponding to center direction
    Ic = findall(x -> abs(imag(x)) > 1e-3 && abs(x) ≈ 1, Me.values)
    # Eigenvalue corresponding to center direction
    Ec = Me.values[Ic[2]];
    # Eigenvector corresponding to center direction
    Vc = Me.vectors[:, Ic[2]];
    if length(Ic) > 2
        Ec = Me.values[Ic[gp["center_idx"]]]
        Vc = Me.vectors[:, Ic[gp["center_idx"]]];
    end

    # INVARIANT CIRCLES
    # Replicate and reshape the vector containing winding angles
    THTj = repeat(thtj', d);
    THTj = THTj[:];
    # Rotation number
    rho = atan(imag(Ec), real(Ec));
    # Correction factor
    cf = exp.((-im*rho*tk)/Tpo);
    # Perturbation from periodic orbit
    uhat = zeros(d*M, N);
    Ulin = zeros(d, M);
    # Construct each invariant circle
    for idx in 1:M
        # Vector index
        vdx = idx*d - (d - 1):idx*d;
        # Time at index
        t = tk[idx];
        # Periodic orbit STM
        phit = sol(t);
        Ulin[:, idx] = phit[1:d];
        phit = reshape(phit[d+1:end], (dp, dp));
        phit = phit[1:d, 1:d];
        # Propagate
        yt1 = cf[idx]*phit*Vc;
        YT1 = repeat(yt1, N, 1)
        # Perturbation
        uh = gp["K"]*(cos.(THTj).*real.(YT1) - sin.(THTj).*imag.(YT1));
        uhat[vdx, :] = reshape(uh, d, N);
    end
    Ulin = repeat(Ulin[:], 1, N);
    Ulin = Ulin + uhat;
    #U = Ulin[:];
    #uhat = uhat[:];

    # ORGANIZE CURVES
    if gp["ctswitch"]
        Uling = zeros(N*d*M, 1);
        uhatg = zeros(N*d*M, 1);
        for i in 1:M
            k = d*N;
            idx = (i*k - (k - 1)):i*k;
            # Invariant curve
            u = Ulin[(i*d - (d - 1)):i*d, :];
            # Store invariant curve
            Uling[idx] = copy(u[:]);
            # Tangent direction
            uh = uhat[(i*d - (d - 1)):i*d, :];
            uhatg[idx] = copy(uh[:]);
        end
        U = Uling;
        uhat = uhatg;
    end

    # INITIAL ANGLE PARTIALS
    DU0t1 = gp["DFr"]*U[1:d*N];
    DU0t0 = zeros(size(DU0t1));
    w = [2*pi/Tpo; rho/Tpo];
    dx = zeros(6);
    for i in 1:gp["N"]
        # Index of torus point
        idx = (i*gp["d"] - (gp["d"] -1)):i*gp["d"];
        # Torus point
        x = U[idx];
        # Update state derivative
        gp["f"](dx, x, p, 0);
        #gp["f"](dx, x, p, (2*pi/gp["N"])*(i - 1));
        # Partial of initial solution wrt tht0
        DU0t0[idx] = (1/w[1])*(dx - w[2]*DU0t1[idx]);
    end

    # INTIAL GUESS CONSTRUCTION
    # Initial guess
    z0 = [U; Tpo; rho; w; p];
    # Initial Tangent Vector
    phi0 = [uhat; zeros(4 + length(p), 1)]/sqrt(dot(uhat, uhat)/N*M);
    # Initial angle partials
    Ud0["DU0t0"] = DU0t0;
    Ud0["DU0t1"] = DU0t1;
    # Initial solution
    Ud0["z0"] = z0;
    # Patch times
    Ud0["pt"] = LinRange(0, Tpo, gp["M"] + 1);
    Ud0["pt"] = 0:step(Ud0["pt"]):Ud0["pt"][end]-1e-5;
    Ud0["Ut"] = zeros(d*N, M);
    # Initialize matrix of STMs
    Ud0["phiT"] = zeros(d*N, d*N, M);
    # Initialize parameters partials 
    Ud0["dxdpt"] = zeros(d*N, dp-d, M)
    # Vector field at final integrated state
    fT = zeros(size([Ud0["Ut"][:, end]; repeat(zeros(dp*dp*N, 1))]));
    Ud0["fT"] = fT[1:d*N];
    
    return (z0, phi0, Ud0)

end

function F_gmos(z, Ud, gp)
    # =================================
    # =================================
    #
    # F - GMOS
    # By: Damennick Henry
    # Date: 11/10/21
    # 
    # Description: Compute the error vector for calculating qpos using the GMOS method
    #
    # Inputs
    #       z - Current solution guess
    #       Ud - Torus solution dictionary
    #       gp - GMOS parameters dictionary
    #
    # Outputs
    #       F - Error vector
    #       Ud - Torus solution dictionary
    #
    # Dependencies:
    #
    # Notes: See  "Fully Numerical Methods for Continuing Families of Quasi-Periodic Invariant Tori in Astrodyanamics" by N. Baresi et al for description of the GMOS algorithm
    #      
    # =================================
    # =================================

    # INITIALIZATION
    # Dictionary extraction
    # State dimension
    d = gp["d"];
    # State and parameter dimension
    dp = gp["dp"]
    # Number of invariant curve points
    N = gp["N"];
    # Number of multiple shooting points
    M = gp["M"];
    # Discrete Fourier transform matrices
    # Fr performs a Fourier transform
    # IFr performs the inverse fourier transfrom
    # DFr computes the partial of a torus function wrt an angle
    Fr = gp["Fr"]; IFr = gp["IFr"]; DFr = gp["DFr"];
    
    # Solution extraction
    # Torus function
    U = z[1:d*N*M];
    # Stroboscopic time
    T = z[d*N*M+1];
    # Rotation number
    rho = z[d*N*M+2]
    # Frequencies
    w0 = z[d*N*M+3];
    w1 = z[d*N*M+4];
    # Parameters
    p = z[d*N*M+5:end];
    
    # INTEGRATION 
    # Initialize matrix of integrated states
    Ud["Ut"] = zeros(d*N, M);
    # Initialize matrix of STMs
    Ud["phiT"] = zeros(d*N, d*N, M);
    # Initialize things
    Ud["dxdpt"] = zeros(d*N, dp-d, M)
    # Integration times
    time = [collect(Ud["pt"]); T];
    # Integrate each patch point
    for i in 1:M
        # Extract invariant circle
        k = d*N;
        idx = (i*k - (k - 1)):i*k;
        u = U[idx];
        # Time span for integration
        ts = (time[i], time[i+1]);
        # Initial STM
        phi0 = 1.0*Matrix(I, dp, dp);
        # Create ODE problem
        prob = ODEProblem(gp["gEOM"], [u; repeat(phi0[:], N)], ts, p);
        # Solve the ODE
        sol = solve(prob, VCABM(), abstol = 1e-16, reltol = 3e-14);
        
        # Final states
        uphit = sol.u[end];
        # Integrated states
        Ud["Ut"][:, i] = uphit[1:d*N];
        # STMS
        phit = uphit[d*N+1:end];
        # Put STMs into a block matrix
        for j in 1:N
            # Block matrix index
            blkidx = (j*d - (d - 1)):j*d;
            # STM index
            stmidx = (j*(dp^2) - ((dp^2) - 1)):j*dp^2;
            # Entire STM
            phiT = reshape(phit[stmidx], dp, dp);
            # State STM
            Ud["phiT"][blkidx, blkidx, i] = phiT[1:d, 1:d];
            # Partial with respect to parameter
            Ud["dxdpt"][blkidx, :, i] = phiT[1:d, end-(dp-d)+1:end];
        end
    end
    # Vector field at final integrated state
    fT = zeros(size([Ud["Ut"][:, end]; repeat(zeros(dp*dp*N, 1))]));
    gp["gEOM"](fT, [Ud["Ut"][:, end]; repeat(zeros(dp*dp*N, 1))], p, T);
    Ud["fT"] = fT[1:d*N];

    # Rotation matrix in the Fourier domain
    rotation_matrix!(gp["Q"], rho, d, N);
    R = IFr*gp["Q"]*Fr;
    
    # Phase conditions
    P = phasecon_gmos(U[1:d*N], Ud, DFr, gp);

    # Frequency consistency constraints
    c0 = T*w0 - 2*pi;
    c1 = T*w1 - rho;

    # Error vector
    F = zeros(d*N*M + length(gp["pconidx"]) + 2);
    # Invariance constraint
    F[1:d*N] = R*Ud["Ut"][:, end] - U[1:d*N];
    # Multiple shooting constraints
    for i in 2:M
        # Invariant circle index
        k = d*N;
        invidx = (i*k - (k-1)):i*k;
        # Patch points
        uptch = U[invidx];
        # Integrated previous points
        F[invidx] = Ud["Ut"][:, i - 1] - uptch;
    end

    # Store phase conditions
    F[d*N*M + 1:d*N*M + length(gp["pconidx"])] .= P[gp["pconidx"]];

    # Store consistency constraints
    F[end-1:end] = [c0; c1];

    # Parameterization constraints
    for i in 1:length(gp["s"])
        si = gp["s"][i];
        si = si(z, Ud);
        append!(F, si);
    end
    
    return (F, Ud)

end

function DF_gmos(z, Ud, gp)
    # =================================
    # =================================
    #
    # DF - GMOS
    # By: Damennick Henry
    # Date: 11/11/21
    # 
    # Description: Compute the error vector Jacobian for calculating qpos using the GMOS method
    #
    # Inputs
    #       z - Current solution guess
    #       Ud - Torus solution dictionary
    #       gp - GMOS parameters dictionary
    #
    # Outputs
    #       F - Error vector
    #       Ud - Torus solution dictionary
    #
    # Dependencies:
    #
    # Notes: See  "Fully Numerical Methods for Continuing Families of Quasi-Periodic Invariant Tori in Astrodyanamics" by N. Baresi et al for description of the GMOS algorithm
    #      
    # =================================
    # =================================

    # INITIALIZATION
    # Dictionary extraction
    # State dimension
    d = gp["d"];
    # State and parameter dimension
    dp = gp["dp"]
    # Number of invariant curve points
    N = gp["N"];
    # Number of multiple shooting points
    M = gp["M"];
    # Fourier matrices
    Fr = gp["Fr"]; IFr = gp["IFr"]; DFr = gp["DFr"];
    
    # Solution extraction
    # Torus function
    U = z[1:d*N*M];
    # Stroboscopic time
    T = z[d*N*M+1];
    # Rotation number
    rho = z[d*N*M+2]
    # Frequencies
    w0 = z[d*N*M+3];
    w1 = z[d*N*M+4];
    # Parameters
    p = z[d*N*M+5:end];

    # Rotation matrix
    rotation_matrix!(gp["Q"], rho, d, N);
    R = IFr*gp["Q"]*Fr;

    # Initialize Jacobian
    DF = zeros(d*N*M + length(gp["pconidx"]) + 2, d*N*M + 4 + length(p))

    # INVARIANCE CONDITION PARTIALS
    # Partial of inv condition wrt U0
    DF[1:d*N, 1:d*N] = -1.0*Matrix(I, d*N, d*N);
    if M == 1
        DF[1:d*N, 1:d*N] += R*Ud["phiT"][:, :, end];
    else
        DF[1:d*N, end-(d*N)-3-length(p):end-4-length(p)] = R*Ud["phiT"][:, :, end];
    end
    # Partial of invariance condition wrt T
    DF[1:d*N, d*N*M + 1] = R*Ud["fT"];
    # Partial of invariance condition wrt rho
    DF[1:d*N, d*N*M + 2] = -DFr*R*Ud["Ut"][:, end];
    # Partial of invariance condition wrt system parameters
    DF[1:d*N, d*N*M + 5:end] = R*Ud["dxdpt"][:, :, end];
    
    # MULTIPLE SHOOTING PARTIALS
    # Loop through patch points
    for i in 2:M
        k = d*N;
        # Row index
        rdx = (i*k - (k - 1)):i*k;
        # Column index
        cdx = rdx .- k;
        # Partials with respect to U0
        DF[rdx, cdx] = Ud["phiT"][:, :, i-1];
        DF[rdx, rdx] = -1.0*Matrix(I, k, k);
        # Partials with respect to system parameters
        DF[rdx, d*N*M + 5:end] = Ud["dxdpt"][:, :, i-1];
    end

    # PHASE CONDITION PARTIALS
    # All of the phase condition partials
    DP = [Ud["DU0t0"]'/N; Ud["DU0t1"]'/N];
    # Grab the appropriate partials
    DF[d*N*M + 1:d*N*M + length(gp["pconidx"]), 1:d*N] = DP[gp["pconidx"], :]

    # CONSISTENCY CONSTRAINT PARTIALS
    # c0
    DF[end-1, d*N*M+1] = w0;
    DF[end-1, d*N*M+3] = T;
    # c1
    DF[end, d*N*M+1] = w1;
    DF[end, d*N*M+2] = -1;
    DF[end, d*N*M+4] = T;

    # PARAMETERIZATION CONSTRAINT PARTIALS
    # Loop through partial functions
    for i in 1:length(gp["Ds"])
        # Partial function i
        Dsi = gp["Ds"][i];
        # Partial function value
        Dsi = Dsi(z, Ud);
        # Append the partial to the Jacobian
        DF = [DF; Dsi];
    end
    

    return (DF, Ud)
end

function Ud_finalization_gmos!(z, Ud, gp)
    # =================================
    # =================================
    #
    # Torus Dictionary Finalization - GMOS
    # By: Damennick Henry
    # Date: 11/11/21
    # 
    # Description: Finalize the torus dictionary at the end of finding a solution for a continuation process
    #
    # Inputs
    #       z - State
    #       Us - Torus solution dictionary
    #
    # Outputs
    #       Ud - Torus solution dictionary
    #
    # Dependencies:
    #
    #      
    # =================================
    # =================================

    # INITIALIZATION
    # Dictionary extraction
    # State dimension
    d = gp["d"];
    # State and parameter dimension
    dp = gp["dp"]
    # Number of invariant curve points
    N = gp["N"];
    # Number of multiple shooting points
    M = gp["M"];
    # Fourier matrices
    Fr = gp["Fr"]; IFr = gp["IFr"]; DFr = gp["DFr"];
    
    # Solution extraction
    # Torus function
    U = @view(z[1:d*N*M]);
    # Stroboscopic time
    T = z[d*N*M+1];
    # Rotation number
    rho = z[d*N*M+2]
    # Frequencies
    w0 = z[d*N*M+3];
    w1 = z[d*N*M+4];
    # Parameters
    p = @view(z[d*N*M+5:end]);

    # Torus function partial wrt theta 1
    mul!(Ud["DU0t1"], DFr, @view(U[1:d*N]));
    # State derivative
    dx = zeros(d);
    for i in 1:N
        # Index of torus point
        idx = (i*d - (d -1)):i*d;
        # Torus point
        x = @view(U[idx]);
        # Partial of initial solution wrt tht0
        gp["f"](dx, x, p, 0);
        #gp["f"](dx, x, p, (2*pi/gp["N"])*(i - 1));
        Ud["DU0t0"][idx] = (1/w0)*(dx - w1*Ud["DU0t1"][idx]);
    end

    # Update previous solutions
    Ud["z0"] = z;
    Ud["w"] = [w0; w1];
    #Ud["pt"] = LinRange(0, T, M + 1);
    #Ud["pt"] = 0:step(Ud["pt"]):Ud["pt"][end]-1e-5;

    return Ud
    
end

function eom_gmos!(xdot, x, p, t, xi, idx, gp)
    # =================================
    # =================================
    #
    # EOM - GMOS
    # By: Damennick Henry
    # Date: 7/13/21
    # 
    # Description: Equations of motion function for computing N points in GMOS
    #
    # Inputs
    #       x - System state
    #       p - System parameters
    #       t - Time
    #       gp - GMOS parameters dictionary
    #
    # Outputs
    #       Xdot - Vector field and variational equations update
    #
    # Dependencies:
    #
    # Notes:
    #      
    # =================================
    # =================================


    # Loop through each state
    
    for i in 1:gp["N"]
        # Set indicies of interest
        idx[i*gp["d"] - (gp["d"] -1):i*gp["d"]] .= 1;
        idx[(i*(gp["dp"]^2) - ((gp["dp"]^2) - 1) + gp["N"]*gp["d"]):((i*(gp["dp"]^2)) + gp["N"]*gp["d"])] .= 1;
        # Get state
        xi[1:gp["d"]] .= @view(x[(i*gp["d"] - (gp["d"] -1)):i*gp["d"]]);
        # Get STM
        xi[gp["d"]+1:end] .= @view(x[(i*(gp["dp"]^2) - ((gp["dp"]^2) - 1) + gp["N"]*gp["d"]):((i*(gp["dp"]^2)) + gp["N"]*gp["d"])]);
        # Get derivative
        gp["eom"](@view(xdot[idx]), xi, p, t);
        #gp["eom"](@view(xdot[idx]), xi, p, t + (2*pi/gp["N"])*(i - 1));
        # Reset indicies
        idx[i*gp["d"] - (gp["d"] -1):i*gp["d"]] .= 0;
        idx[(i*(gp["dp"]^2) - ((gp["dp"]^2) - 1) + gp["N"]*gp["d"]):((i*(gp["dp"]^2)) + gp["N"]*gp["d"])] .= 0;
    end

    return nothing
end

function compute_full_torus!(z, Ud, gp)
    # =================================
    # =================================
    #
    # Compute Full Torus
    # By: Damennick Henry
    # Date: 11/16/21
    # 
    # Description: Compute the entire torus from an invariant curve. All of the tangent directions are also computed
    #
    # Inputs
    #       z - Torus solution vector
    #       Ud - Torus solution dictionary
    #       gp - GMOS parameters dictionary
    #
    # Outputs
    #       Ud - Torus solution dictionary updated with the full torus information      
    # 
    #
    # Dependencies:
    #
    # Notes: 
    #      
    # =================================
    # =================================

    # INITIALIZATION
    # State dimension
    d = gp["d"];
    # State and parameter dimension
    dp = gp["dp"];
    # Number of invariant curve points
    N = gp["N"];
    # Number of multiple shooting points
    M = gp["M"];
    # Fourier matrices
    Fr = gp["Fr"]; IFr = gp["IFr"]; DFr = gp["DFr"];
    # Desired angles in the tht0 direction
    tht0 = gp["tht0"];
    # Stable and unstable tangent directions for the invariant curve
    vs = Ud["vs"]; vu = Ud["vu"];
    # Torus function
    U = z[1:d*N*M];
    # Stroboscopic time
    T = z[d*N*M+1];
    # Frequencies
    w0 = z[d*N*M+3];
    w1 = z[d*N*M+4];
    # Parameters
    p = z[d*N*M+5:end];

    # INTEGRATE INVARIANT CURVE
    # Initial STM
    phi0 = 1.0*Matrix(I, dp, dp);
    # Create ODE Problem
    prob = ODEProblem(gp["gEOM"], [U[1:N*d]; repeat(phi0[:], N)], (0.0, T), p);
    # Integrate the curve
    sol = solve(prob, VCABM(), abstol = 1e-16, reltol = 3e-14);

    # ROTATE POINTS AND TANGENT DIRECTIONS
    # Initialize torus function matrix and full tangent direction matrices
    U = zeros(d*N, length(tht0));
    Vs = zeros(d*N, length(tht0)); Vu = zeros(d*N, length(tht0));
    # Loop through each desired angle
    for i in 1:length(tht0)
        # Time corresponding to tht0
        t = tht0[i]/w0;
        # Rotation of tht1
        tht1 = t*w1;
        # Rotation Matrix
        rotation_matrix!(gp["Q"], tht1, d, N);
        # Rotation transformation
        R = IFr*gp["Q"]*Fr;
        # Invariant circle at time T
        u = sol(t);
        U[:, i] = R*u[1:d*N];
        phit = u[d*N+1:end];
        
        for j in 1:N
            # STM index
            stmidx = j*(dp^2) - ((dp^2) - 1):j*(dp^2);
            PHIT = reshape(phit[stmidx], dp, dp);
            PHIT = PHIT[1:d, 1:d];
            vst = PHIT*vs[(j*d - (d-1)):j*d];
            vut = PHIT*vu[(j*d - (d-1)):j*d];
            vst = vst/norm(vst);
            vut = vut/norm(vut);
            Vs[(j*d - (d-1)):j*d, i] = copy(vst);
            Vu[(j*d - (d-1)):j*d, i] = copy(vut);
        end
        Vs[:, i] = R*Vs[:, i];
        Vu[:, i] = R*Vu[:, i];
    end

    # Store information
    Ud["U"] = U; Ud["Vs"] = Vs; Ud["Vu"] = Vu;
end

function phasecon_gmos(U, Ud, DFr, gp)
    # =================================
    # =================================
    #
    # Phase condition - GMOS
    # By: Damennick Henry
    # Date: 11/11/21
    # 
    # Description: Compute the phase condition for calculating qpos using the GMOS method
    #
    # Inputs
    #       U - Matrix of invariant circle points
    #       Ud - Torus solution dictionary
    #       DFr - Fourier partial computation matrix
    #       gp - GMOS parameters dictionary
    #
    # Outputs
    #       P - Vector of phase conditions
    #
    # Dependencies:
    #
    # Notes: See  "Fully Numerical Methods for Continuing Families of Quasi-Periodic Invariant Tori in Astrodyanamics" by N. Baresi et al for description of the GMOS algorithm
    #      
    # =================================
    # =================================

    # INITIALIZATION
    # State dimension
    d = gp["d"];
    # Number of invariant curve points
    N = gp["N"];
    # Previous solution
    z0 = Ud["z0"];
    # Previous partials with respect to the torus angles
    DU0t0 = Ud["DU0t0"];
    DU0t1 = Ud["DU0t1"];

    # PHASE CONDITIONS
    p0 = (1/N)*(U - z0[1:d*N])'*DU0t0;
    p1 = (1/N)*(U - z0[1:d*N])'*DU0t1;
    # Schilder phase conditions
    #p0 = (1/N)*U'*DU0t0; 
    #p1 = (1/N)*U'*DU0t1;
    P = [p0; p1];

    return P
end
