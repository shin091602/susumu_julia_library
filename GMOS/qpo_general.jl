# =================================
# =================================
#
# QPO - General
# By: Damennick Henry
# Date: 6/25/21
# 
# Description: Define the functions used for all QPO computation algorithms
#
# Dependencies:
#      
# =================================
# =================================

mutable struct qpo_solution
    # State dimension
    d
    # Vector field
    f
    # Vector field Jacobian
    dfdx
    # Equations of motion function
    EOM
    # Number of grid points
    N
    # Torus function grid
    U
    # Previous solution
    z0
    # Previous solution partials
    DU0t0
    DU0t1
    # Solution partials
    DUt0
    DUt1
    # State at time T
    UT
    # State transiton matrices at time T
    phiT
    # Vector field at time T
    fT
    # Patch times for multiple shooting
    pt
    # Stable tangent direction
    Vs
    # Unstable tangent
    Vu
    # Dummy "data struct" for continuation
    data
end

Base.copy(Us::qpo_solution) = qpo_solution(Us.d, Us.f, Us.dfdx, Us.EOM, Us.N, Us.U, Us.z0, Us.DU0t0, Us.DU0t1, Us.DUt0, Us.DUt1, Us.UT, Us.phiT, Us.fT, Us.pt, Us.Vs, Us.Vu, Us.data)

function initqpoguess_po!(Us, zpo, DR, K)
    # =================================
    # =================================
    #
    # Initial QPO Guess - Periodic Orbit
    # By: Damennick Henry
    # Date: 6/25/21
    # 
    # Description: Generate an initial QPO guess from the center tangent direction of a periodic orbit
    #
    # Dependencies:
    #
    # Notes: See  "Fully Numerical Methods for Continuing Families of Quasi-Periodic Invariant Tori in Astrodyanamics" by N. Baresi et al
    #      
    # =================================
    # =================================

    # INITIALIZATION 
    # State dimension
    d = Us.d
    # Torus grid
    N = Us.N;
    # State and period of the periodic orbit
    xpo = zpo[1:6];
    Tpo = zpo[7];
    # Time vector
    tk = collect(LinRange(0, Tpo, Us.N[1] + 1)); pop!(tk);
    # Winding angle
    thtj = collect(LinRange(0, 2*pi, Us.N[2] + 1)); pop!(thtj);

    # COMPUTE MONODROMY MATRIX
    # Initial STM
    phi0 = 1.0*Matrix(I, d, d);
    # Initial condition
    X0 = [xpo; phi0[:]];
    ts = (0, Tpo);
    # Create ODE Problem
    prob = ODEProblem(Us.EOM, X0, ts);
    # Solve ODE
    sol = solve(prob, VCABM(), abstol = 1e-16, reltol = 3e-14);

    # Final state solution
    xT = sol.u[end];
    # Monodromy matrix
    phiT = reshape(xT[d + 1:end], (d, d));

    # INVARIANT CIRCLE
    # Eigendata of monodromy matrix
    Me = eigen(phiT);
    # Eigenvalue index corresponding to center direction
    Ic = findall(x -> abs(imag(x)) > 1e-3 && abs(x) ≈ 1, Me.values)
    # Eigenvalue corresponding to center direction
    Ec = Me.values[Ic[2]];
    # Eigenvector corresponding to center direction
    Vc = Me.vectors[:, Ic[2]];
    # INITIAL GUESS
    # Replicate and reshape the vector containing winding angles
    THTj = repeat(thtj', d);
    THTj = THTj[:];
    # Rotation number
    rho = atan(imag(Ec), real(Ec));
    # Correction factor
    cf = exp.((-im*rho*tk)/Tpo);
    # Perturbation from periodic orbit
    uhat = zeros(d*N[1], N[2]);
    Ulin = zeros(d, N[1]);
    for idx in 1:N[1]
        # Vector index
        vdx = idx*d - (d - 1):idx*d;
        # Time at index
        t = tk[idx];
        # Periodic orbit STM
        phit = sol(t);
        Ulin[:, idx] = phit[1:d];
        phit = reshape(phit[d+1:end], (d, d));
        # Propagate
        yt1 = cf[idx]*phit*Vc;
        YT1 = repeat(yt1, N[2], 1)
        # Perturbation
        uh = K*(cos.(THTj).*real.(YT1) - sin.(THTj).*imag.(YT1));
        uhat[vdx, :] = reshape(uh, d, N[2]);
    end

    Ulin = repeat(Ulin[:], 1, N[2]);
    Ulin = Ulin + uhat;
    
    Us.z0 = [Ulin[:]; (2*pi)/Tpo; rho/Tpo];
    DU = DU_dft(Ulin, DR, Us.d);
    Us.DU0t0 = DU[1];
    Us.DU0t1 = DU[2];

    #xdx = 1:6:length(Ulin[:]);
    #ydx = 2:6:length(Ulin[:]);
    #zdx = 3:6:length(Ulin[:]);
    #display(scatter3d(Us.z0[xdx], Us.z0[ydx], Us.z0[zdx]))
    return (Us, uhat)
    
end

function fourier_matrix(d, N)
    # =================================
    # =================================
    #
    # Fourier Matrix
    # By: Damennick Henry
    # Date: 7/6/21
    # 
    # Description: Generate a Fourier matrices for an invariant circle
    #
    # Inputs
    #       d - State dimension
    #       N - Number of points on the invariant circle
    #
    # Outputs 
    #       R - Fourier transformation matrix
    #       IR - Inverse Fourier transformation matrix
    #       DR - Partial with respect to the angle transformation matrix
    #
    # Dependencies:
    #
    # Notes: Copied over from MATLAB code written by Nicola Baresi
    #      
    # =================================
    # =================================

    # Initialize
    D = zeros(d*N, d*N) .+ 0im;
    R = zeros(d*N, d*N) .+ 0im;
    IR = zeros(d*N, d*N);
    DR = zeros(d*N, d*N);

    # Order
    K = (-(N-1)/2):1:((N-1)/2);
    J = findall(k -> k == 0, K)
    J = J[1];
    K = [K[J:end]; K[1:J-1]];
    
    for i in 1:N
        # k index
        k = 1;
        # Row index
        rdx = d*(i-1)+1:d*i;
        for j in 1:N
            # Column index
            cdx = d*(j-1)+1:d*j;
            # DFT matrix
            D[rdx, cdx] = 1.0*Matrix(I, d, d).*exp(2*pi*K[i]*K[j]*1.0*im/N);
            # Inverse Fourier matrix
            if mod(j, 2) == 0
                k += 1;
                IR[rdx, cdx] = 1.0*cos(2*pi*(i - 1)*K[k]/N)*Matrix(I, d, d);
                DR[rdx, cdx] = 1.0*-K[k]*sin(2*pi*(i - 1)*K[k]/N)*Matrix(I, d, d);
            else
                IR[rdx, cdx] = 1.0*sin(2*pi*(i - 1)*K[k]/N)*Matrix(I, d, d);
                DR[rdx, cdx] = 1.0*K[k]*cos(2*pi*(i - 1)*K[k]/N)*Matrix(I, d, d);
            end
        end
        # A0
        IR[rdx, 1:d] = 1.0*Matrix(I, d, d);
    end

    # Complex -> Real transformation
    k = 0
    for i in 2:2:(N - 1)
        # Column index
        k += 1;
        # Indexing
        r1 = (d*(i - 1) + 1):d*i;
        r2 = (d*i + 1):d*(i + 1);
        c1 = (d*k + 1):d*(k + 1);
        c2 = (d*(N - k) + 1):d*(N - k + 1);
        # Define R
        R[r1, c1] = 1.0*Matrix(I, d, d)/N;
        R[r1, c2] = 1.0*Matrix(I, d, d)/N;
        R[r2, c1] = -1.0*im*Matrix(I, d, d)/N;
        R[r2, c2] = 1.0*im*Matrix(I, d, d)/N;
    end
    R[1:d, 1:d] = 1.0*Matrix(I, d, d)/N;
    R = R*D;
    if abs(maximum(imag(R))) > 1e-16
        @warn "R is not exactly real"
    end
    R = real(R);

    # Derivative matrix 
    DR = DR*R;
    if abs(maximum(imag(DR))) > 1e-16
        @warn "DR is not exactly real"
    end
    DR = real(DR);

    return (R, IR, DR)
    #=
    file = matopen("/home/dam/Documents/04_CSML/03_Code/01_MATLAB/07_HeteroclinicConnections/03_EML1toL2/jDR.mat", "w");
    write(file, "foo", DR)
    close(file)
    =#
    
end

function rotation_matrix(rho, d, N)
    # =================================
    # =================================
    #
    # Rotation Matrix
    # By: Damennick Henry
    # Date: 7/6/21
    # 
    # Description: Generate a Fourier matrices for an invariant circle
    #
    # Inputs
    #       rho - Negative of angle to rotate by
    #       d - State dimension
    #       N - Number of points on the invariant circle
    #
    # Outputs 
    #       Q - Rotation matrix in Fourier domain
    #
    #
    # Dependencies:
    #
    # Notes: Copied over from MATLAB code written by Nicola Baresi
    #      
    # =================================
    # =================================

    # Indicies
    K = (-(N-1)/2):1:((N-1)/2);
    J = findall(k -> k == 0, K)
    J = J[1];
    K = [K[J:end]; K[1:J-1]];

    # Initialize
    Q = zeros(d*N, d*N);
    Id = 1.0*Matrix(I, d, d)
    k = 1;
    for i in 2:2:(N -1)
        # Column index
        k += 1;
        # Indexing
        row1 = (d*(i - 1) + 1):d*i
        row2 = d*i+1:d*(i + 1);
        # Rotating Fourier coefficients
        mul!(@view(Q[row1, row1]), Id, cos(K[k]*rho));
        mul!(@view(Q[row1, row2]), Id, -sin(K[k]*rho));
        mul!(@view(Q[row2, row1]), Id, sin(K[k]*rho));
        mul!(@view(Q[row2, row2]), Id, cos(K[k]*rho));
    end
    # Add first block
    Q[1:d, 1:d] = Id;
    return Q
end

function rotation_matrix!(Q, rho, d, N)
    # =================================
    # =================================
    #
    # Rotation Matrix
    # By: Damennick Henry
    # Date: 7/6/21
    # 
    # Description: Generate a Fourier matrices for an invariant circle
    #
    # Inputs
    #       rho - Negative of angle to rotate by
    #       d - State dimension
    #       N - Number of points on the invariant circle
    #
    # Outputs 
    #       Q - Rotation matrix in Fourier domain
    #
    #
    # Dependencies:
    #
    # Notes: Copied over from MATLAB code written by Nicola Baresi
    #      
    # =================================
    # =================================

    # Indicies
    K = (-(N-1)/2):1:((N-1)/2);
    J = findall(k -> k == 0, K)
    J = J[1];
    K = [K[J:end]; K[1:J-1]];

    Id = 1.0*Matrix(I, d, d);
    k = 1;
    for i in 2:2:(N -1)
        # Column index
        k += 1;
        # Indexing
        row1 = (d*(i - 1) + 1):d*i
        row2 = d*i+1:d*(i + 1);
        # Rotating Fourier coefficients
        mul!(@view(Q[row1, row1]), Id, cos(K[k]*rho));
        mul!(@view(Q[row1, row2]), Id, -sin(K[k]*rho));
        mul!(@view(Q[row2, row1]), Id, sin(K[k]*rho));
        mul!(@view(Q[row2, row2]), Id, cos(K[k]*rho));
    end
    # Add first block
    Q[1:d, 1:d] = Id;
    return nothing
end

function DU_dft(U, par)
    # =================================
    # =================================
    #
    # DU - DFT
    # By: Damennick Henry
    # Date: 7/6/21
    # 
    # Description: Compute parital derivatives
    #
    # Inputs
    #       U - Matrix of torus points
    #       DR - Partial derivative 
    #
    # Outputs 
    #       DU - Derivative 
    #
    # Dependencies:
    #
    # Notes: 
    #      
    # =================================
    # =================================

    # INITIALIZATION
    # State dimension
    d = par["d"]
    # Number of torus points in each dimension
    N = par["NT"];
    # Fourier matrices
    DFr0 = par["DFr0"]; DFr1 = par["DFr1"];
    # Partial with respect to theta 1
    DUt1 = DFr1*U;
    DUt0 = copy(DUt1);
    for i in 1:N[2]
        # Vector index
        vdx = (d*(i - 1) + 1):d*i;
        # Invariant circle wrt theta 0
        ivc0 = U[vdx, :];
        # Partial with respect to theta 0
        DUt0[vdx, :] = reshape(DFr0*ivc0[:], d, N[1]);
    end
    
    return (DUt0, DUt1)

end

function DU_dft!(DUt0, DUt1, U, par)
    # =================================
    # =================================
    #
    # DU - DFT
    # By: Damennick Henry
    # Date: 7/6/21
    # 
    # Description: Compute parital derivatives
    #
    # Inputs
    #       U - Matrix of torus points
    #       DR - Partial derivative 
    #
    # Outputs 
    #       DU - Derivative 
    #
    # Dependencies:
    #
    # Notes: 
    #      
    # =================================
    # =================================

    # INITIALIZATION
    # State dimension
    d = par["d"]
    # Number of torus points in each dimension
    N = par["NT"];
    # Fourier matrices
    DFr0 = par["DFr0"]; DFr1 = par["DFr1"];
    # Partial with respect to theta 1
    mul!(DUt1, DFr1, U);
    # Partial wrt theta 0
    for i in 1:N[2]
        # Vector index
        vdx = (d*(i - 1) + 1):d*i;
        # Invariant circle wrt theta 0
        ivc0 = @view(U[vdx, :]);
        # Partial with respect to theta 0
        DUt0[vdx, :] = reshape(DFr0*@view(ivc0[:]), d, N[1]);
    end
    
    return nothing

end

function dWdU_dft(W, p)
    # =================================
    # =================================
    #
    # dWdU - DFT
    # By: Damennick Henry
    # Date: 7/6/21
    # 
    # Description: Compute parital derivatives of whisker wrt angles
    #
    # Inputs
    #       U - Matrix of torus points
    #       DR - Partial derivative 
    #
    # Outputs 
    #       DU - Derivative 
    #
    # Dependencies:
    #
    # Notes: 
    #      
    # =================================
    # =================================
    dWdt0 = copy(W);
    dWdt1 = copy(W);
    (~, npts) = size(W)
    for i in 1:npts
        (DUt0, DUt1) = DU_dft(reshape(W[:, i], p["NT"][2]*p["d"], p["NT"][1]), p);
        dWdt0[:, i] = DUt0[:]; dWdt1[:, i] = DUt1[:];
    end

    return (dWdt0, dWdt1)
end

function interpolate_torus(U, N, d)
    # =================================
    # =================================
    #
    # Interpolate Torus
    # By: Damennick Henry
    # Date: 7/10/21
    # 
    # Description: Interpolate a grid of torus function points
    #
    # Inputs
    #       U - Matrix of torus points
    #       N - Number of torus points on each grid
    # Outputs 
    #       F - Array of interpolants
    #
    # Dependencies:
    #
    # Notes:
    #      
    # =================================
    # =================================
    # Collect grid points
    tht0 = LinRange(-2*pi, 4*pi, N[1]*3 + 1);
    tht0 = tht0[1]:step(tht0):tht0[end]-1e-5;
    tht1 = LinRange(-2*pi, 4*pi, N[2]*3 + 1);
    tht1 = tht1[1]:step(tht1):tht1[end]-1e-5;
    THTs = (tht0, tht1)
    # Interpolate the torus
    F = fill([], 1, d)
    for dim in 1:d
        # Index for dimension
        i = dim:d:N[2]*d;
        # Values for dimension
        vi = repeat(U[i, :]', 3, 3);
        # Interpolant for dimension
        Fi = CubicSplineInterpolation(THTs, vi, extrapolation_bc = NaN);
        # Interpolant
        F[dim] = [Fi];

    end
    return F
end

function interpolate_whisker(W, t, N, d)
    # =================================
    # =================================
    #
    # Interpolate Whisker
    # By: Damennick Henry
    # Date: 8/19/21
    # 
    # Description: Interpolate a grid of torus function points
    #
    # Inputs
    #       W - Matrix of whisker points
    #       t - Time at each point
    #       N - Number of torus points on each grid
    #       d - State dimension
    # Outputs 
    #       F - Array of interpolants
    #
    # Dependencies:
    #
    # Notes: 
    #      
    # =================================
    # =================================

    # Collect grid points
    tht0 = LinRange(-2*pi, 4*pi, N[1]*3 + 1);
    tht0 = tht0[1]:step(tht0):tht0[end]-1e-5;
    tht1 = LinRange(-2*pi, 4*pi, N[2]*3 + 1);
    tht1 = tht1[1]:step(tht1):tht1[end]-1e-5;
    i = sortperm(t);
    t = t[i];
    W = W[:, i];
    t = t[1]:t[2]-t[1]:t[end];
    THTs = (tht0, tht1, t)
    nt = length(t);
    
    # Interpolate the torus
    F = Array{Interpolations.FilledExtrapolation, 1}(undef, d);
    vi = zeros(N[1]*3, N[2]*3, nt)
    for dim in 1:d
        # Index for dimension
        i = dim:d:N[1]*N[2]*d;
        # Values for dimension
        for j in 1:nt
            vi[:, :, j] = repeat(reshape(W[i, j], N[1], N[2])', 3, 3);
        end
        # Interpolant for dimension
        F[dim] = CubicSplineInterpolation(THTs, vi, extrapolation_bc = NaN);
        #Fi = interpolate((collect(THTs[1]), collect(THTs[2]), collect(THTs[3])), vi, Gridded(Linear()));
        #println(typeof(Fi))
        # Interpolant
        #F[dim] = [Fi];
    end

    return F
end

function plot_torus(F, d)
    # =================================
    # =================================
    #
    # Interpolate Torus
    # By: Damennick Henry
    # Date: 7/10/21
    # 
    # Description: Interpolate a grid of torus function points
    #
    # Inputs
    #       F - Array of interpolants
    # Outputs 
    #
    # Dependencies:
    #
    # Notes: 
    #      
    # =================================
    # =================================

    # Angles to plot
    tht0 = LinRange(0, 2*pi, 100);
    tht1 = LinRange(0, 2*pi, 100);

    # Position vectors
    r = zeros(length(tht0), length(tht1), round(Int, d/2))
    for i in 1:round(Int, d/2)
        Fi = F[i];
        Fi = Fi[1];
        r[:, :, i] = [Fi(t0, t1) for t0 in tht0, t1 in tht1];
    end
    
    surf(r[:, :, 1], r[:, :, 2], r[:, :, 3])
end

function rotate_torus(U, tht, par)
    # =================================
    # =================================
    #
    # Rotate Torus
    # By: Damennick Henry
    # Date: 8/17/21
    # 
    # Description: Rotate torus grid U by tht
    #
    # Inputs
    #       U - Torus grid
    #       tht - Angle to rotate over
    #       par - Parameters dictionary
    #
    # Outputs 
    #       Ur - Rotated torus grid
    #
    # Dependencies:
    #
    # Notes: 
    #      
    # =================================
    # =================================

    # INITIALIZE 
    # Parameters extraction
    # State dimension
    d = par["d"]
    # Number of torus points in each dimension
    N = par["NT"];
    # Fourier matrices
    Fr0 = par["Fr0"]; IFr0 = par["IFr0"]; DFr0 = par["DFr0"];
    Fr1 = par["Fr1"]; IFr1 = par["IFr1"]; DFr1 = par["DFr1"];

    # ROTATE FIRST ANGLE 
    # Rotation matrix
    Q = rotation_matrix(-tht[1], d, N[1]);
    R = IFr0*Q*Fr0;
    # Reorder the grid
    Ut = zeros(N[1]*d, N[2]);
    for j in 1:N[2]
        vjdx = (j*d - (d - 1)):j*d;
        for i in 1:N[1]
            vidx = (i*d - (d - 1)):i*d;
            Ut[vidx, j] = U[vjdx, i];
        end
    end
    # Rotate
    Ur = R*Ut;
    Ut = zeros(N[2]*d, N[1]);
    # Reorder the grid
    for i in 1:N[1]
        vidx = (i*d - (d - 1)):i*d;
        for j in 1:N[2]
            vjdx = (j*d - (d - 1)):j*d;
            Ut[vjdx, i] = Ur[vidx, j];
        end
    end
    Ur = Ut;

    # ROTATE SECOND ANGLE
    # Rotation matrix
    Q = rotation_matrix(-tht[2], d, N[2]);
    R = IFr1*Q*Fr1;
    Ur = R*Ur;
    
    return Ur
end

function wmatch(z, Us, wd)

    return z[end-1] - wd
end

function Dwmatch(z, Us)
    # Initialize
    Jwd = zeros(1, length(z));
    # Jacobian
    Jwd[end-1] = 1;
    return Jwd
end

function intmatch(z, int_func, N, d, Ct)
    Csd = 0;
    for i in 1:N
        C = int_func(z[(i*d - (d-1)):i*d])
        Csd += C;
    end

    Csd = (1/prod(N))*Csd - Ct;
    return Csd;
end

function Dintmatch(z, Dint_func, N, d)
    DCsd = zeros(1, length(z));
    for i in 1:N
        vdx = (i*d - (d-1)):i*d;
        DCsd[vdx] = Dint_func(z[vdx]);
    end
    DCsd = (1/prod(N))*DCsd;
end