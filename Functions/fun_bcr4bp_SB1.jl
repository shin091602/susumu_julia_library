function fun_bcr4bp_SB1!(dx, x, parameter, t)
    """
    This function computes the derivatives of position and velocity for an object under the influence of two primary bodies 
    in the Sun-B1 rotating frame.
    """
    # Parameter unpacking
    mu_SB1 = parameter[1]
    mu_EM = parameter[2]
    a_S = parameter[3]
    a_EM = parameter[4]
    w_M = parameter[5]
    theta_M0 = parameter[6]

    # Sun - B1 rotating frame
    # t       : non-dimensional time
    # x       : non-dimensional position and velocity, x = [x y z vx vy vz]'
    # mu_SB1   : mass ratio in the Sun - B1 frame, muSB1 = (m1 + m2) / (m1 + m2 + ms)
    # mu_EM    : mass ratio in the Earth - Moon frame, muEM = m2 / (m1 + m2)
    # a_S      : non-dimensional Sun orbit radius
    # a_EM     : non-dimensional distance between Earth and Moon
    # w_M      : non-dimensional Moon angular velocity
    # theta_M0 : initial value of the Moon's angle

    # Calculate the Moon's angle (initial angle + angular velocity * time)
    theta = theta_M0 + w_M * t

    # The distances to the primary bodies
    r_s = sqrt((x[1] + mu_SB1)^2 + x[2]^2 + x[3]^2)
    # Earth position
    r_E = [-mu_EM * cos(theta), -mu_EM * sin(theta), 0] .* (a_EM / a_S) .+ [1 - mu_SB1, 0, 0]
    # Moon position
    r_M = [(1 - mu_EM) * cos(theta), (1 - mu_EM) * sin(theta), 0] .* (a_EM / a_S) .+ [1 - mu_SB1, 0, 0]
    r1 = norm([x[1] - r_E[1], x[2] - r_E[2], x[3] - r_E[3]])
    r2 = norm([x[1] - r_M[1], x[2] - r_M[2], x[3] - r_M[3]])
    r1 = sqrt((x[1] - r_E[1])^2 + (x[2] - r_E[2])^2 + (x[3] - r_E[3])^2)
    r2 = sqrt((x[1] - r_M[1])^2 + (x[2] - r_M[2])^2 + (x[3] - r_M[3])^2)

    # Derivatives of position and velocity
    dx[1] = x[4]
    dx[2] = x[5]
    dx[3] = x[6]
    dx[4] =  2 * x[5]   + x[1]  - (1 - mu_SB1) * (x[1] + mu_SB1) / r_s^3    - mu_SB1 * (1 - mu_EM) * (x[1] - r_E[1]) / r1^3     - mu_SB1 * mu_EM * (x[1] - r_M[1]) / r2^3
    dx[5] = -2 * x[4]   + x[2]  - (1 - mu_SB1) * x[2] / r_s^3               - mu_SB1 * (1 - mu_EM) * (x[2] - r_E[2]) / r1^3     - mu_SB1 * mu_EM * (x[2] - r_M[2]) / r2^3
    dx[6] =                     - (1 - mu_SB1) * x[3] / r_s^3               - mu_SB1 * (1 - mu_EM) * (x[3] - r_E[3]) / r1^3     - mu_SB1 * mu_EM * (x[3] - r_M[3]) / r2^3

end
