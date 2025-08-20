function fun_bcr4bp_EMS!(dx, x, parameter, t)
    """
    Computes the derivatives of position and velocity in the Earth-Moon rotating frame 
    with the influence of the Sun.
    """

    # Parameter unpacking
    mu = parameter[1]
    m_S = parameter[2]
    a_S = parameter[3]
    omega_S = parameter[4]
    theta_S0 = parameter[5] 

    # x  : non-dimensional position and velocity, x = [x y z vx vy vz]'
    # mu : mass ratio of the primaries 
    # m_S : Nondimensional mass of the Sun
    # a_S : Nondimensional Sun orbit radius
    # omega_S : Nondimensional Sun angular velocity
    # theta_S0 : Initial Sun angle
    # t  : non-dimensional true anomaly

    # Calculate the Sun's angle (initial angle + angular velocity * time)
    theta_S = theta_S0 + omega_S * t

    # Compute the Sun's X and Y coordinates on its circular orbit
    Xs = a_S * cos(theta_S)
    Ys = a_S * sin(theta_S)
    Zs = 0.0 # The Sun does not move in the Z direction 

    # The distances to the primary bodies
    # Distance from the spacecraft to the Earth
    r1 = sqrt((x[1] + mu)^2 + x[2]^2 + x[3]^2)
    # Distance from the spacecraft to the Moon
    r2 = sqrt((x[1] - 1 + mu)^2 + x[2]^2 + x[3]^2)
    # Distance from the spacecraft to the Sun
    r3 = sqrt((x[1] - Xs)^2 + (x[2] - Ys)^2 + (x[3] - Zs)^2)

    # Derivatives of position and velocity
    dx[1] = x[4]
    dx[2] = x[5]
    dx[3] = x[6]
    dx[4] =  2*x[5] +x[1]  -(1-mu)*(x[1]+mu)/r1^3 -mu*(x[1]-1+mu)/r2^3 -m_S*(x[1]-Xs)/r3^3 -m_S/a_S^3*Xs
    dx[5] = -2*x[4] +x[2]  -(1-mu)*x[2]/r1^3      -mu*x[2]/r2^3        -m_S*(x[2]-Ys)/r3^3 -m_S/a_S^3*Ys
    dx[6] =                -(1-mu)*x[3]/r1^3      -mu*x[3]/r2^3        -m_S*(x[3]-Zs)/r3^3 -m_S/a_S^3*Zs

end
