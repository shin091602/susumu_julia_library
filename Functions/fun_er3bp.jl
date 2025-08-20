function fun_er3bp!(dx, x, parameter, t)
    """
    This function calculates the derivatives of position and velocity for an object under the influence of two primary bodies
    in the ER3BP framework. The results are stored directly in the `dx` array to avoid memory allocation.
    """
    
    # Parameter unpacking 
    mu = parameter[1]
    e = parameter[2]

    # x  : non-dimensional position and velocity, x = [x y z vx vy vz]'
    # mu : mass ratio of the primaries (parameter[1])
    # e  : eccentricity (parameter[2])
    # t  : non-dimensional true anomaly

    # The distances to the primary bodies
    r1 = sqrt((mu + x[1])^2 + x[2]^2 + x[3]^2)
    r2 = sqrt((1 - mu - x[1])^2 + x[2]^2 + x[3]^2)

    # Partial derivatives of the effective potential
    dUdx = x[1] - (1 - mu) * (x[1] + mu) / r1^3 - mu * (x[1] - 1 + mu) / r2^3
    dUdy = x[2] - (1 - mu) * x[2] / r1^3 - mu * x[2] / r2^3
    dUdz = x[3] - (1 - mu) * x[3] / r1^3 - mu * x[3] / r2^3

    # Define the derivatives of position and velocit
    dx[1] = x[4]
    dx[2] = x[5]
    dx[3] = x[6]
    dx[4] = 2 * x[5] + dUdx / (1 + e * cos(t))
    dx[5] = -2 * x[4] + dUdy / (1 + e * cos(t))
    dx[6] = -x[3] + dUdz / (1 + e * cos(t))
end
