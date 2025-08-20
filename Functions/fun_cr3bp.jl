function fun_cr3bp!(dx, x, parameter, t) # !は入力変数を書き換えてしまう（例えばdx）破壊的な関数につける
    """
    This function calculates the derivatives of position and velocity for an object under the influence of two primary bodies
    in the CR3BP framework. The results are stored directly in the `dx` array to avoid memory allocation.
    """
    # Parameter unpacking 
    mu = parameter[1]

    # t   : non-dimensional time
    # x   : non-dimensional position and velocity, x = [x y z vx vy vz]'
    # mu  : mass ratio of the primaries 

    # Calculate the distances to the primary bodies
    r1 = sqrt((mu + x[1])^2 + x[2]^2 + x[3]^2)
    r2 = sqrt((x[1] - 1 + mu)^2 + x[2]^2 + x[3]^2)

    # Define the derivatives of position and velocity
    dx[1] = x[4]
    dx[2] = x[5]
    dx[3] = x[6]
    dx[4] =  2*x[5] +x[1] -(1-mu)*(x[1]+mu)/r1^3 -mu*(x[1]-1+mu)/r2^3
    dx[5] = -2*x[4] +x[2] -(1-mu)*x[2]/r1^3      -mu*x[2]/r2^3
    dx[6] =               -(1-mu)*x[3]/r1^3      -mu*x[3]/r2^3
end
