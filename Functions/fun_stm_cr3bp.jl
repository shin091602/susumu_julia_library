function fun_stm_cr3bp!(dx, x, parameter, t)

    mu = parameter[1]

    # Equations of motion
    r1 = sqrt((mu + x[1])^2 + x[2]^2 + x[3]^2)
    r2 = sqrt((1 - mu - x[1])^2 + x[2]^2 + x[3]^2)
    
    # Define the derivatives of position and velocity
    dx[1] = x[4]
    dx[2] = x[5]
    dx[3] = x[6]
    dx[4] =  2*x[5] +x[1] -(1-mu)*(x[1]+mu)/r1^3 -mu*(x[1]-1+mu)/r2^3
    dx[5] = -2*x[4] +x[2] -(1-mu)*x[2]/r1^3      -mu*x[2]/r2^3
    dx[6] =               -(1-mu)*x[3]/r1^3      -mu*x[3]/r2^3
    
    # STM
    A = fun_A_cr3bp(x[1:6], mu)
    phi = reshape(x[7:end], 6, 6)
    dphi = A * phi
    dx[7:42] = dphi[:]
    
end
