# =================================
# =================================
#
# Periodic Orbits - Multiple Shooting
# By: Damennick Henry
# Date: 2/1/22
# 
# Description: Define functions for computing peroidic orbits using multiple shooting
#
# Dependencies:
#
# Notes:
#
# =================================
# =================================

function F_poms(z, Xd, msp)
    # =================================
    # =================================
    #
    # F - Periodic Orbit Single Shooting
    # By: Damennick Henry
    # Date: 2/1/22
    # 
    # Description: Compute error vector for calculating a periodic orbit using single shooting.
    #
    # Inputs:
    #       z - Free variable vector [x; T; p]
    #       Xd - Periodic orbit solution dictionary
    #       msp - Shooting parameters dict
    #
    # Outputs:
    #       F  - Error vector
    #       Xd - Periodic orbit solution dictionary
    #
    # Dependencies:
    #
    # Notes:
    #
    # =================================
    # =================================

    # INITIALIZATION
    # State dimension
    d = msp["d"];
    # Number of patches
    M = msp["M"];
    # Patch points
    x = z[1:d*M];
    # Period
    T = z[d*M+1];
    # Parameters
    p = z[d*M+2:end];
    # Parameter dimension
    dp = length(p);

    # INTEGRATION
    # Initial STM
    phi0 = 1.0*Matrix(I, d + dp, d + dp);
    # Patch times
    pt = LinRange(0, T, M + 1); pt = 0:step(pt):pt[end]-1e-5;
    Xd["xT"] = zeros(d, M);
    Xd["phiT"] = zeros(d+ dp, d + dp, M)
    # "Full state"
    Xd["Xfull"] = [];
    # Integrate each patch point
    for i = 1:M
        # Initial condition
        X0 = [x[(i-1)*d + 1:i*d]; phi0[:]];
        # Time span for integration
        ts = (pt[i], i == M ? T : pt[i+1]);
        # Create ODE Problem
        prob = ODEProblem(msp["eom"], X0, ts, p);
        # Solve ODE
        sol = solve(prob, VCABM(), abstol = 1e-16, reltol = 3e-14);
        # Final state
        Xd["xT"][:, i] = sol.u[end][1:d];
        # Final STM
        Xd["phiT"][:, :, i] = reshape(sol.u[end][d+1:end], d + dp, d + dp);
        # Vector field at final state
        msp["f"](@view(Xd["fT"][:, i]), Xd["xT"][:, i], p, i == M ? T : pt[i+1]);
        # "Full state"
        if Xd["Xfull"] != []
            Xd["Xfull"] = [Xd["Xfull"] sol[1:d, 2:end]];
            Xd["tfull"] = [Xd["tfull"]; sol.t[2:end]]
        else
            Xd["Xfull"] = sol[1:d, :];
            Xd["tfull"] = sol.t;
        end
    end

    #display(plot(Xd["Xfull"][1, :], Xd["Xfull"][2, :]))

    F = zeros(d*M)
    for i = 2:M
        j = i - 1;
        jdx = (j-1)*d + 1:j*d;
        F[jdx] = Xd["xT"][:, j] - x[(i-1)*d + 1:i*d];
    end

    F[(M-1)*d + 1:M*d] = Xd["xT"][:, end] - x[1:d];
    

    if msp["phsflg"]
        # Phase condition
        phs = dot(z[1:6] - Xd["x0"], Xd["f0"]);
        # Error vector
        F = [F; phs];
    end

    # Parameterization constraints
    for i in 1:length(msp["s"])
        si = msp["s"][i];
        si = si(z, Xd);
        append!(F, si);
    end
    return (F, Xd)
end

function DF_poms(z, Xd, msp)
    # =================================
    # =================================
    #
    # DF - Periodic Orbit Multiple Shooting
    # By: Damennick Henry
    # Date: 2/1/22
    # 
    # Description: Compute error vector Jacobian for calculating a periodic orbit using multiple shooting.
    #
    # Inputs:
    #       z - Free variable vector [x; T; p]
    #       Xd - Periodic orbit solution dictionary
    #       msp - Shooting parameters dict
    #
    # Outputs:
    #       F  - Error vector
    #       Xd - Periodic orbit solution dictionary
    #
    # Dependencies:
    #
    # Notes:
    #
    # =================================
    # =================================

    # INITIALIZATION
    # State dimension
    d = msp["d"];
    # Number of patches
    M = msp["M"];
    # Patch points
    x = z[1:d*M];
    # Period
    T = z[d*M+1];
    # Parameters
    p = z[d*M+2:end];
    # Parameter dimension
    dp = length(p);

    DF = zeros(d*M, d*M + 1 + dp)

    for c = 2:M
        # Row index
        r = c - 1;
        rdx = (r-1)*d + 1:r*d;
        # Column index
        cdx = (c-1)*d + 1:c*d;
        

        # Partial wrt integrated patch point
        DF[rdx, rdx] = Xd["phiT"][1:d, 1:d, r];
        # Partial wrt patch point
        DF[rdx, cdx] = -1.0*Matrix(I, d, d);
        # Partial wrt period
        DF[rdx, d*M + 1] = Xd["fT"][:, r]/M;
        # Partial wrt parameters 
        DF[rdx, d*M + 2:end] = Xd["phiT"][1:d, d+1:end, r];
    end

    # PERIODICITY PARTIAL
    # Partial wrt integrated patch point
    DF[(M-1)*d + 1:M*d, (M-1)*d + 1:M*d] = Xd["phiT"][1:d, 1:d, end];
    # Partial wrt patch point
    DF[(M-1)*d + 1:M*d, 1:d] = -1.0*Matrix(I, d, d);
    # Partial wrt period
    DF[(M-1)*d + 1:M*d, d*M + 1] = Xd["fT"][:, end]/M;
    # Partial wrt parameters 
    DF[(M-1)*d + 1:M*d, d*M + 2:end] = Xd["phiT"][1:d, d+1:end, end];

    # PHASE CONSTRAINT PARTIAL
    if msp["phsflg"]
        # Phase condition
        Dphs = zeros(1, d*M + 1 + dp);
        Dphs[1, 1:6] = Xd["f0"]';
        # Error vector
        DF = [DF; Dphs];
    end

    # Loop through partial functions
    for i in 1:length(msp["Ds"])
        # Partial function i
        Dsi = msp["Ds"][i];
        # Partial function value
        Dsi = Dsi(z, Xd);
        # Append the partial to the Jacobian
        DF = [DF; Dsi];
    end


    return (DF, Xd)
end

function Xd_finalization_ms!(z, Xd, msp)
    # =================================
    # =================================
    #
    # Xd Finalization - Multiple Shooting
    # By: Damennick Henry
    # Date: 11/14/21
    # 
    # Description: Finalize the multiple shooting dictionary
    #
    # Inputs:
    #       z - Free variable vector [x; T; p]
    #       Xd - Periodic orbit solution dict
    #
    # Outputs:
    #       DF - Error vector Jacobian
    #       Xd - Periodic orbit state dict
    #
    # Dependencies:
    #
    # Notes: See Baresi disseration for details on algorithm.
    #
    # =================================
    # =================================

    Xd["x0"] = z[1:msp["d"]];
    msp["f"](Xd["f0"], Xd["x0"], z[msp["d"]+2:end], 0);

    return Xd
end

function po_ms_initialization(z0, p)
    # =================================
    # =================================
    #
    # Periodic Orbit Multiple Shooting Initialization
    # By: Damennick Henry
    # Date: 3/7/22
    # 
    # Description: Initialize an initial PO guess for multiple shooting
    #
    # Inputs:
    #       z0 - Free variable vector [x; T; p]
    #       Xd0 - Periodic orbit solution dict
    #
    # Outputs:
    #       
    #
    # Dependencies:
    #
    # Notes: See Baresi disseration for details on algorithm.
    #
    # =================================
    # =================================

    # Patch times
    pt = LinRange(0, z0[7], p["M"] + 1); pt = 0:step(pt):pt[end]-1e-5;
    # Integrate initial guess
    prob = ODEProblem(p["f"], z0[1:6], (0.0, z0[7]), z0[8:end]);
    # Solve ODE
    sol = solve(prob, VCABM(), abstol = 1e-16, reltol = 3e-14,saveat = pt);
    #display(plot([u[1] for u in sol.u], [u[2] for u in sol.u], [u[3] for u in sol.u]))
    #adsf
    # Collect patch points
    x = sol[:, :];
    # Create initial guess
    z0 = [x[:]; z0[7]; z0[8:end]]

    # Update dict
    Xd0 = Dict();
    Xd0["fT"] = zeros(p["d"], p["M"])
    Xd0["x0"] = z0[1:p["d"]];
    Xd0["f0"] = zeros(p["d"]); p["f"](Xd0["f0"], Xd0["x0"], 0, 0);

    return (z0, Xd0)
end

function po_refine(z0, Xd0, p)
    # =================================
    # =================================
    #
    # Periodic Orbit Refine
    # By: Damennick Henry
    # Date: 3/7/22
    # 
    # Description: Refine a PO solution to a certain value
    #
    # Inputs:
    #       z0 - Free variable vector [x; T; p]
    #       Xd0 - Periodic orbit solution dict
    #
    # Outputs:
    #       z0 - Refined solution
    #
    # Dependencies:
    #
    # Notes: See Baresi disseration for details on algorithm.
    #
    # =================================
    # =================================

    # Refine the initial guess
    # Define error vector, error vector Jacobian and dictionary finalizaiton function
    F(z, Xd) = F_poms(z, Xd, p);
    DF(z, Xd) = DF_poms(z, Xd, p);
    vpf(z, Xd) = Xd_finalization_ms!(z, Xd, p);
    # Continuation parameters
    pac_p = pac_parameters(0, 1e-12, 20, 0, 0, 0, true, 0, 0, false, true, false);
    display(z0)
    # Peform pseudo-arclength continuation
    Zpo = pa_cont(z0, 1e-5, Xd0, [], F, DF, vpf, pac_p);
    if Zpo != []
        return Zpo.z
    else
        return []
    end
end

