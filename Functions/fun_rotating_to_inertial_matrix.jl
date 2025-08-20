function fun_rotating_to_inertial_matrix(C, dtheta_dt)
    """
    This function constructs a transformation matrix from a rotating frame to an inertial frame.
    """
    # C  : Rotation Matrix
    # dtheta_dt  : Angular Velocity

    rotating_matrix = [
        # Transformation of the position vector (top 3 rows)
        C[1,1] C[1,2] C[1,3] 0 0 0;
        C[2,1] C[2,2] C[2,3] 0 0 0;
        C[3,1] C[3,2] C[3,3] 0 0 0;
        
        # Transformation of the velocity vector (bottom 3 rows)
        dtheta_dt * C[1,2] -dtheta_dt * C[1,1] 0 C[1,1] C[1,2] C[1,3];
        dtheta_dt * C[2,2] -dtheta_dt * C[2,1] 0 C[2,1] C[2,2] C[2,3];
        dtheta_dt * C[3,2] -dtheta_dt * C[3,1] 0 C[3,1] C[3,2] C[3,3]
    ]
    return rotating_matrix
end
