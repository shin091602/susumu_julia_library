function fun_cr3bp_parameter(y)
    """
    This function returns four parameters (mass ratio, primary-secondary distance, secondary radius, and angular velocity) 
    for a specified system in the circular restricted three-body problem (CR3BP).
    """

    if y == 1
        # Sun-Earth
        # Mass ratio
        mu = 3.00245E-06
        # Primary-secondary distance [km] 
        a_1 = 1.49598E+08 
        # Secondary radius [-]
        a_s = 4.25874E-05 
        # Angular velocity [rad/s]
        w_1 = 1.99097E-07 

    elseif y == 2
        # Earth-Moon
        mu = 1.21536E-02
        a_1 = 3.84400E+05
        a_s = 4.52003E-03
        w_1 = 2.66167E-06

    elseif y == 3
        # Sun-Mars
        mu = 3.22605E-07
        a_1 = 2.27943824E+08
        a_s = 1.48699E-05
        w_1 = 1.0585760E-07

    elseif y == 4
        # Mars-Phobos
        mu = 1.66100E-08
        a_1 = 9.37600E+03
        a_s = 1.18387E-03
        w_1 = 2.28040E-04

    elseif y == 5
        # Mars-Deimos
        mu = 2.30046E-09
        a_1 = 2.34580E+04
        a_s = 2.34580E+04
        w_1 = 2.34580E+04

    elseif y == 6
        # Sun-Jupiter
        mu = 9.5479E-04
        a_1 = 7.7841E+08
        a_s = 9.1844E-04
        w_1 = 1.6800E-08

    else
        error("Invalid body ID. Use a value from 1 to 6.")
    end

    return mu, a_1, a_s, w_1
end
