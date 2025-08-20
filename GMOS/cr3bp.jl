# =================================
# =================================
#
# CR3BP
# By: Damennick Henry
# Date: 6/23/21
# 
# Description: Define various functions for the CR3BP
#
# Dependencies:
#
# Notes: 
#
# =================================
# =================================


function f_CR3BP!(f, x, mu)
    # =================================
    # =================================
    #
    # f - CR3BP
    # By: Damennick Henry
    # Date: 6/23/21
    # 
    # Description: Compute the vector field for the circular restricted three body problem
    #
    # Inputs:
    #       x - Spacecraft state
    #       mu - CR3BP mass fraction
    #
    # Outputs:
    #       f - Vector field at x
    #
    # Dependencies:
    #
    # Notes: 
    #
    # =================================
    # =================================
    
    # Velocitites
    f[1] = x[4]; 
    f[2] = x[5]; 
    f[3] = x[6]; 

    # Accelerations
    f[4] = x[1]+(-1)*mu*((-1)+mu+x[1])*(((-1)+mu+x[1])^2+x[2]^2+x[3]^2)^(-3/2)+((-1)+mu)*(mu+x[1])*((mu+x[1])^2+x[2]^2+x[3]^2)^(-3/2)+2*x[5]; 
    f[5] = x[2]+(-1)*mu*x[2]*(((-1)+mu+x[1])^2+x[2]^2+x[3]^2)^(-3/2)+((-1)+mu)*x[2]*((mu+x[1])^2+x[2]^2+x[3]^2)^(-3/2)+(-2)*x[4]; 
    f[6] = (-1)*mu*x[3]*(((-1)+mu+x[1])^2+x[2]^2+x[3]^2)^(-3/2)+((-1)+mu)*x[3]*((mu+x[1])^2+x[2]^2+x[3]^2)^(-3/2); 

    return nothing
end

function Df_CR3BP!(dfdx, x, mu, d)
    # =================================
    # =================================
    #
    # Df - CR3BP
    # By: Damennick Henry
    # Date: 6/23/21
    # 
    # Description: Compute the vector field Jacobian for the circular restricted three body problem
    #
    # Inputs:
    #       x - Spacecraft state
    #       mu - CR3BP mass fraction
    #
    # Outputs:
    #       dfdx - Vector field Jacobian at x
    #
    # Dependencies:
    #
    # Notes: 
    #
    # =================================
    # =================================

    dfdx[1:d, 1:d] .= @SMatrix [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; 1+3*mu*((-1)+mu+x[1])^2*((( -1)+mu+x[1])^2+x[2]^2+x[3]^2)^(-5/2)+(-1)*mu*(((-1)+mu+x[1])^2+x[2]^2+x[3]^2)^(-3/2)+3*(1+(-1)*mu)*(mu+x[1])^2*((mu+x[1])^2+x[2]^2+x[3]^2)^(-5/2)+(-1)*(1+(-1)*mu)*((mu+x[1])^2+x[2]^2+x[3]^2)^(-3/2) 3*mu*((-1)+mu+x[1])*x[2]*(((-1)+mu+x[1])^2+x[2]^2+x[3]^2)^(-5/2)+3*(1+(-1)*mu)*(mu+x[1])*x[2]*((mu+x[1])^2+x[2]^2+x[3]^2)^(-5/2) 3*mu*((-1)+mu+x[1])*x[3]*(((-1)+mu+x[1])^2+x[2]^2+x[3]^2)^(-5/2)+3*(1+(-1)*mu)*(mu+x[1])*x[3]*((mu+x[1])^2+x[2]^2+x[3]^2)^(-5/2) 0 2 0; 3*mu*((-1)+mu+x[1])*x[2]*(((-1)+mu+x[1])^2+x[2]^2+x[3]^2)^(-5/2)+3*(1+(-1)*mu)*(mu+x[1])*x[2]*((mu+x[1])^2+x[2]^2+x[3]^2)^(-5/2) 1+3*mu*x[2]^2*(((-1)+mu+x[1])^2+x[2]^2+x[3]^2)^(-5/2)+(-1)*mu*(((-1)+mu+x[1])^2+x[2]^2+x[3]^2)^(-3/2)+3*(1+(-1)*mu)*x[2]^2*((mu+x[1])^2+x[2]^2+x[3]^2)^(-5/2)+(-1)*(1+(-1)*mu)*((mu+x[1])^2+x[2]^2+x[3]^2)^(-3/2) 3*mu*x[2]*x[3]*(((-1)+mu+x[1])^2+x[2]^2+x[3]^2)^(-5/2)+3*(1+(-1)*mu)*x[2]*x[3]*((mu+x[1])^2+x[2]^2+x[3]^2)^(-5/2) (-2) 0 0; 3*mu*((-1)+mu+x[1])*x[3]*(((-1)+mu+x[1])^2+x[2]^2+x[3]^2)^(-5/2)+3*(1+(-1)*mu)*(mu+x[1])*x[3]*((mu+x[1])^2+x[2]^2+x[3]^2)^(-5/2) 3*mu*x[2]*x[3]*(((-1)+mu+x[1])^2+x[2]^2+x[3]^2)^(-5/2)+3*(1+(-1)*mu)*x[2]*x[3]*((mu+x[1])^2+x[2]^2+x[3]^2)^(-5/2) 3*mu*x[3]^2*(((-1)+mu+x[1])^2+x[2]^2+x[3]^2)^(-5/2)+(-1)*mu*(((-1)+mu+x[1])^2+x[2]^2+x[3]^2)^(-3/2)+3*(1+(-1)*mu)*x[3]^2*((mu+x[1])^2+x[2]^2+x[3]^2)^(-5/2)+(-1)*(1+(-1)*mu)*((mu+x[1])^2+x[2]^2+x[3]^2)^(-3/2) 0 0 0]; 

    return nothing
end

function jacobi_integral(x, mu)
    # =================================
    # =================================
    #
    # Jacobi Integral
    # By: Damennick Henry
    # Date: 6/24/21
    # 
    # Description: Compute the Jacobi integral for the circular restricted three body problem
    #
    # Inputs:
    #       x - Spacecraft state
    #       mu - CR3BP mass fraction
    #
    # Outputs:
    #       C - Vector field at x
    #
    # Dependencies:
    #       LinearAlgebra
    # Notes: 
    #
    # =================================
    # =================================

    r = x[1:3];
    v = x[4:6];

    U = (1 - mu)/norm(r + mu*[1; 0; 0]) + mu/norm(r + (mu-1)*[1; 0; 0]) + 0.5*(r[1]^2 + r[2]^2);    

    C = 2*U - norm(v)^2;

    return C
end

function jacobi_integral(x, mu)
    # =================================
    # =================================
    #
    # Jacobi Integral
    # By: Damennick Henry
    # Date: 6/24/21
    # 
    # Description: Compute the Jacobi integral for the circular restricted three body problem
    #
    # Inputs:
    #       x - Spacecraft state
    #       mu - CR3BP mass fraction
    #
    # Outputs:
    #       C - Jacobi integral at x
    #
    # Dependencies:
    #
    # Notes: 
    #
    # =================================
    # =================================

    r = x[1:3];
    v = x[4:6];

    U = (1 - mu)/norm(r + mu*[1; 0; 0]) + mu/norm(r + (mu-1)*[1; 0; 0]) + 0.5*(r[1]^2 + r[2]^2);    

    C = 2*U - norm(v)^2;

    return C
end

function DC(x, mu)
    # =================================
    # =================================
    #
    # DC
    # By: Damennick Henry
    # Date: 6/24/21
    # 
    # Description: Compute the gradient of the Jacobi integral for the circular restricted three body problem
    #
    # Inputs:
    #       x - Spacecraft state
    #       mu - CR3BP mass fraction
    #
    # Outputs:
    #       dCdx - Jacobi integral gradient at x
    #
    # Dependencies:
    #
    # Notes: 
    #
    # =================================
    # =================================

    dC = [2*(x[1]+(-1)*mu*((-1)+mu+x[1])*(((-1)+mu+x[1])^2+x[2]^2+x[3]^2)^(-3/2)+((-1)+mu)*(mu+x[1])*((mu+x[1])^2+x[2]^2+x[3]^2)^(-3/2)),2*(x[2]+(-1)*mu*x[2]*(((-1)+mu+x[1])^2+x[2]^2+x[3]^2)^(-3/2)+((-1)+mu)*x[2]*((mu+x[1])^2+x[2]^2+x[3]^2)^(-3/2)),(-2)*mu*x[3]*(((-1)+mu+x[1])^2+x[2]^2+x[3]^2)^(-3/2)+2*((-1)+mu)*x[3]*((mu+x[1])^2+x[2]^2+x[3]^2)^(-3/2),(-2)*x[4],(-2)*x[5],(-2)*x[6]];

    return dC
end

function DU(r, mu)
    # =================================
    # =================================
    #
    # DU
    # By: Damennick Henry
    # Date: 6/25/21
    # 
    # Description: Compute the gradient of the gravitational potential in the planar CR3BP problem
    #
    # Inputs:
    #       r - Two dimensional position vector
    #       mu - CR3BP mass fraction
    #
    # Outputs:
    #       dUdx - Gradient of the gravitational potential
    #
    # Dependencies:
    #
    # Notes: 
    #
    # =================================
    # =================================

    x = r[1];
    y = r[2];

    dUdx = x - (mu*(-1 + mu + x))/(((-1 + mu + x)^2 + y^2)^(3/2))- ((1 - mu)*(mu + x))/(((mu + x)^2 + y^2)^(3/2));

    return dUdx
end

function librationpoints(mu)
    # =================================
    # =================================
    #
    # Equilibrium points
    # By: Damennick Henry
    # Date: 6/25/21
    # 
    # Description: Compute the equilibrium points in the CR3BP
    #
    # Inputs:
    #       mu - CR3BP mass fraction
    #
    # Outputs:
    #       Lix - Libration point x coordinates
    #
    # Dependencies:
    #       DU
    #       Roots
    # Notes: 
    #
    # =================================
    # =================================

    # Tolerance
    tol = 1e-12;

    f(x) = DU([x; 0], mu);
    L3x = find_zero(f, (-20, -mu - tol), Bisection());
    L1x = find_zero(f, (-mu + tol, 1 - mu - tol), Bisection());
    L2x = find_zero(f, (1 - mu + tol, 20), Bisection())

    return (L1x, L2x, L3x)
end

