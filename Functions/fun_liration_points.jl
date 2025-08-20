using Polynomials

function fun_libration_points(mu)
    """
    This function calculates the equilibrium point coordinates in the circular restricted three-body problem (CR3BP). 
    The input variable `mu` represents the mass parameter, and the function outputs the positions of the five Lagrange points: L1, L2, L3, L4, and L5.
    """

    # Define a polynomial (Be careful that the order of coefficients in Julia is the reverse of MATLAB's `roots` function.)
    l = 1 - mu
    p_L1 = Polynomial([ mu^3-l^3, mu^2*l^2+2*( l^2+mu^2), 2*mu*l*(l-mu)+mu-l, l^2-4*l*mu+mu^2, 2*(mu-l), 1.0])
    p_L2 = Polynomial([-mu^3-l^3, mu^2*l^2+2*( l^2-mu^2), 2*mu*l*(l-mu)-mu-l, l^2-4*l*mu+mu^2, 2*(mu-l), 1.0])
    p_L3 = Polynomial([ mu^3+l^3, mu^2*l^2+2*(-l^2+mu^2), 2*mu*l*(l-mu)+mu+l, l^2-4*l*mu+mu^2, 2*(mu-l), 1.0])

    # L1: Find the roots of a quintic equation
    L1roots = Polynomials.roots(p_L1)
    L1_x = only(filter(x -> imag(x) == 0, L1roots) .|> real)
    L1 = [L1_x, 0.0, 0.0]

    # L2: Find the roots of a quintic equation
    L2roots = Polynomials.roots(p_L2)
    L2_x = only(filter(x -> imag(x) == 0, L2roots) .|> real)
    L2 = [L2_x, 0.0, 0.0]

    # L3: Find the roots of a quintic equation
    L3roots = Polynomials.roots(p_L3)
    L3_x = only(filter(x -> imag(x) == 0, L3roots) .|> real)
    L3 = [L3_x, 0.0, 0.0]

    # L4: Calculate algebraically
    L4 = [0.5-mu,  0.5*sqrt(3), 0.0]

    # L5: Calculate algebraically
    L5 = [0.5-mu, -0.5*sqrt(3), 0.0]
    
    return L1, L2, L3, L4, L5
end
