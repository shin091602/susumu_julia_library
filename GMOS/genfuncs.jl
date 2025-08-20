# =================================
# =================================
#
# General Functions 
# By: Damennick Henry
# Date: 7/1/21
# 
# Description: General functions I use
#
# Dependencies:
#
# Notes: See Baresi disseration for details on algorithm.
#
# =================================
# =================================

function eom!(xdot, x, p, t, A, f!::Function, Df!::Function, d)
    # =================================
    # =================================
    #
    # Equations of Motion 
    # By: Damennick Henry
    # Date: 6/25/21
    # 
    # Description:  The equations of motion and variational equations
    #
    # Inputs:
    #       x - State
    #       t - Time
    #       f - Vector field function
    #       Df - Vector field Jacobian function
    #       d - State dimension
    # Outputs:
    #       xdot - Derivative of state/variational equations
    #
    # Dependencies:
    #
    # Notes: See Baresi disseration for details on algorithm.
    #
    # =================================
    # =================================

    # Vector field Jacobian
    Df!(A, x, p, t);
    # State and paramter dimension
    dp = maximum(size(A))
    # STM
    phi = reshape(@view(x[d+1:end]), dp, dp);
    dphi = reshape(@view(xdot[d+1:end]), dp, dp);
    # STM update
    mul!(@view(dphi[1:d, 1:dp]), A, phi);

    # State update
    f!(xdot, x, p, t);
    xdot[d+1:end] .= @view(dphi[:]);
    return nothing
end

function DF_fd(z, vp, F)
    (Fz, vp) = F(z, vp);

    DF = zeros((length(Fz), length(z)));
    dz = 1e-8;
    for idx in 1:length(z)
        dZ = zeros(size(z));
        dZ[idx] = dz;
        (Fdz, vpc) = F(z + dZ, vp);
        DF[:, idx] = (Fdz - Fz)/dz;
    end

    return DF
end

function DF_fd!(z, vp, F!, Fv)
    println("=================================")
    println("ENTERING FINITE DIFFERENCING LOOP")
    println("=================================")
    vpc = deepcopy(vp);
    F!(Fv, z, vpc);
    Fz = copy(Fv);
    DF = zeros((length(Fz), length(z)));
    dz = 1e-8;
    for idx in 1:length(z)
        if idx == 6
            dz = 1e-12;
        end
        dZ = zeros(size(z));
        dZ[idx] = dz;
        vpc = deepcopy(vp);
        zn = z + dZ;
        
        println(string("z comp: ", idx))
        F!(Fv, z + dZ, vpc);
        
        DF[:, idx] = (Fv - Fz)/dz;
    end

    return DF
end


function weighted_normalization!(v, w)

    nd = ndims(v)
    if nd > 1
    (n, m) = size(v)
    else
        m = 1
    end
    for i in 1:m
        v[:, i] = v[:, i]/sqrt(dot(v[:, i].*w, v[:, i]))
    end
end

function normcol!(x)
    for col in eachcol(x)
        col ./= norm(col)
    end
    return x
end