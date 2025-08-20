# =================================
# =================================
#
# Psuedo-arclength Continuation
# By: Damennick Henry
# Date: 6/28/21
# 
# Description: Define functions for performing pseudo-arclength continuation
#
# Dependencies:
#
# Notes:
#
# =================================
# =================================

struct pac_parameters
    # Number of solutions
    n
    # Error Tolerance
    tol
    # Max iterations
    itmax
    # Optimal number of iterations
    optit
    # Maximum step length
    smax
    # Jacobian rows to consider for inital family tangent
    jcr
    # Sparse Jacobian flag
    issparse
    # Number of continuation function points
    N
    # Continuation function point indices
    fpidx
    # Plot flag
    plotflg
    # Correction only flag
    coflg
    # Finite difference check
    fdcheck
end

mutable struct pac_solution
    # Solution vector
    z
    # Varying parameters
    vp
    # Previous solution 
    z0
    # Tangent vector
    s
    # Tangent basis
    phi
end

function pa_cont(z0, s0, vp0, phi0, F, DF, vpf, p)
    # =================================
    # =================================
    #
    # Psuedo-arclength Continuation
    # By: Damennick Henry
    # Date: 6/28/21
    # 
    # Description: Perform pseudo-arclength continuation
    #
    # Inputs
    #       z0 - Initial Guess
    #       s0 - Initial step length size
    #       phi0 - Tangent vector
    #       vp0 - Varying Parameters
    #       F - Error vector function
    #       DF - Error vector Jacobian function
    #       vpf - Varying parameters finalization function
    #       p - Parameters struct      
    #       
    # Outputs
    #
    # Dependencies:
    #
    # Notes:
    #
    # =================================
    # =================================

    # Compute error vector and jacobian
    (Fv, vp0) = F(z0, vp0);
    (DFv, vp0) = DF(z0, vp0);
    # Initial family tangent
    if phi0 == [] && !p.coflg
        phi0 = nullspace(DFv[p.jcr, :]);
        phi0 = phi0/sqrt((dot(phi0[1:p.fpidx], phi0[1:p.fpidx])/p.N) + dot(phi0[p.fpidx+1:end], phi0[p.fpidx+1:end]));
        if phi0[end] < 0
            println("SWITCHING PHI")
            phi0 *= -1;
        end
    end
    z = pac_solution(z0, vp0, z0, s0, phi0);
    # Correct the initial guess
    z = pa_corrector!(z, p, F, DF, vpf);
    if !p.coflg
        println("Corrected")
    end
    #=
    (DFv, vp) = DF(z.z, z.vp);
    svdDF = svd(DFv, full=true);
    z.phi = svdDF.Vt[end, :];
    w = ones(length(z.z));
    w[1:p.fpidx] = w[1:p.fpidx]/p.N;
    weighted_normalization!(z.phi, w)
    
    if length(phi0)  == length(z.phi)
        if acos(dot(z.phi, phi0)/(norm(z.phi)*norm(phi0))) > π/2
            z.phi = -z.phi;
        end
    end
    =#
    
    if p.coflg
        return z
    end
    
    #=svdDF = svd(DFv);
    z.phi = svdDF.Vt[end, :];
    normalize!(z.phi)
    if acos(dot(z.phi, phi0)) > π/2
        z.phi = -z.phi;
    end
    =#
    (DFv, vp) = DF(z.z, z.vp);
    z.phi = nullspace(DFv, atol=1e-12);
    if length(phi0) == length(z.phi)
        if acos(dot(z.phi, phi0)/(norm(z.phi)*norm(phi0))) > π/2
            z.phi = -z.phi;
        end
    end
    #z.phi = -z.phi
    #=
    if z.phi[end] < 0
        println("SWITCHING PHI")
        z.phi = -z.phi;
    end
    =#
    # Store solutions
    Z = Array{pac_solution, 1}(undef, p.n);
    Z[1] = pac_solution(z.z, z.vp, z.z0, z.s, z.phi);
    
    for idx in 2:p.n
        # Predict the next solution
        z.z0 = copy(z.z);
        z.z = z.z0 + z.phi*z.s;
        # Correct the solution
        z = pa_corrector!(z, p, F, DF, vpf);
        #println(idx)
        if z.z == []
            Z = Z[1:idx-1];
            break
        end
        z.vp = copy(z.vp);
        Z[idx] = pac_solution(z.z, deepcopy(z.vp), z.z0, z.s, z.phi);
        Z[idx].vp = deepcopy(z.vp);
        if p.plotflg
            
            #U = z.z[1:25*6];
            #xdx = 1:6:length(U);
            #ydx = 2:6:length(U);
            #zdx = 3:6:length(U);
            if idx == 2
                display(plot(z.vp["Xfull"][1, :], z.vp["Xfull"][2, :], z.vp["Xfull"][3, :], legend = false, aspect_ratio = 1))
                #display(scatter(U[xdx], U[ydx], legend = false, aspect_ratio = 1))
            elseif idx < p.n -1
                display(plot!(z.vp["Xfull"][1, :], z.vp["Xfull"][2, :], z.vp["Xfull"][3, :], legend = false,  aspect_ratio = 1))
                #display(scatter!(U[xdx], U[ydx], legend = false, aspect_ratio = 1))
            else
                display( plot!(z.vp["Xfull"][1, :], z.vp["Xfull"][2, :], z.vp["Xfull"][3, :], legend = false,  aspect_ratio = 1))
            end
        end
    end

    return Z
end

function pa_corrector!(z::pac_solution, p::pac_parameters, F, DF, vpf)
    # =================================
    # =================================
    #
    # Psuedo-arclength Corrector
    # By: Damennick Henry
    # Date: 6/28/21
    # 
    # Description: Correct to a solution
    #
    # Inputs
    #       z0 - Initial Guess
    #       vp - Varying Parameters
    #       phi - Tangent vector
    #       s - pseudo-arclength step size
    #       F - Error vector function
    #       DF - Error vector Jacobian function
    #       vpf - Varying parameter finalization function
    #       p - Parameters struct
    #          
    # Outputs
    #       z - Solution
    #       
    #
    # Outputs
    #
    # Dependencies:
    #
    # Notes:
    #
    # =================================
    # =================================

    # Initial error
    (Fval, z.vp) = F(z.z, z.vp);
    #Fval = [Fval; dot(z.phi, z.z - z.z0) - z.s];
    if !p.coflg
        paccon = (1/p.N)*dot(z.z[1:p.fpidx] - z.z0[1:p.fpidx], z.phi[1:p.fpidx]) + dot(z.z[p.fpidx+1:end]-z.z0[p.fpidx+1:end], z.phi[p.fpidx+1:end]) - z.s;
        Fval = [Fval; paccon];
    end
    
    err = norm(Fval);
    # Iteration counter
    it = 0;

    # Correct guess
    while err > p.tol && it < p.itmax
        # Compute Jacobian
        (DFval, z.vp) = DF(z.z, z.vp);

        if p.fdcheck
            DFfd = DF_fd(z.z, z.vp, F)
            file = matwrite("jDFs.mat", Dict(
                "jDF" => DFval,
                "jDFfd" => DFfd
            ));
            alkds2fw
            return
        end
        
        
        
        # Add pseudo-arclength Jacobian
        if !p.coflg
            Dpaccon = [z.phi[1:p.fpidx]'/p.N z.phi[p.fpidx+1:end]']
            DFval = [DFval; Dpaccon];
        end
        # Correction
        if p.issparse
            dz = sparse(-DFval)\Fval;
        else
            dz = -DFval\Fval;
        end
        # Apply Correction
        z.z = z.z + dz;
        # Error
        (Fval, z.vp) = F(z.z, z.vp);
        #Fval = [Fval; dot(z.phi, z.z - z.z0) - z.s];
        if !p.coflg
            paccon = (1/p.N)*dot(z.z[1:p.fpidx] - z.z0[1:p.fpidx], z.phi[1:p.fpidx]) + dot(z.z[p.fpidx+1:end]-z.z0[p.fpidx+1:end], z.phi[p.fpidx+1:end]) - z.s;
            Fval = [Fval; paccon];
        end
        
        err = norm(Fval);
        #println(err)
        # Update iteration counter
        it = it + 1;
    end

    if err > p.tol
        println("Solution not found")
        z.z = []
        return z
    else
        if p.coflg == false
            println("PAC solution converged")
            println(string("Number of iterations: ", it))
            println(string("Final error:", err))
            println(string("Step size: ", z.s))
            println(z.z[end])
        end
        if !p.coflg
            (DFval, z.vp) = DF(z.z, z.vp);
            Dpaccon = [z.phi[1:p.fpidx]'/p.N z.phi[p.fpidx+1:end]']
            DFval = [DFval; Dpaccon];
            phi = pac_null(DFval, p);
        
            #if dot(phi, z.phi) < 0
            #    println("Switching phi")
            #    phi = -phi;
            #end
            z.phi = phi;
            # Update step length
            Eps = abs(p.optit/it);
            if Eps > 2 || it == 0
                Eps = 2;
            elseif Eps < 0.5
                Eps = 0.5
            end
            if z.s > 0
                z.s = min(p.smax, z.s*Eps);
            else
                z.s = max(p.smax, z.s*Eps);
            end
        end
        z.vp = vpf(z.z, z.vp);
        return z
    end
end

function pac_null(DFv, p)
    # Size of Jacobian
    (n, m) = size(DFv);
    # Tangent basis
    if p.issparse
        phi = sparse(DFv)\[zeros(n - 1, 1); 1];
    else
        phi = DFv\[zeros(n - 1, 1); 1];
    end
    phi = phi/sqrt((dot(phi[1:p.fpidx], phi[1:p.fpidx])/p.N) + dot(phi[p.fpidx+1:end], phi[p.fpidx+1:end]));
    return phi
end