function condition(x, t, integrator)
    """
    This function is used to identify the moment when the trajectory of the manifold crosses 
    the plane perpendicular to the x-axis at x = 1 - μ, which passes through the secondary primary body.
    """
    # Retrieving the value of mu
    mu = integrator.p

    # Conditions for evaluation  
    # mode =  1 : All cases where x = 1 - μ and dx > 0 are satisfied  
    # mode =  0 : All cases where x = 1 - μ is satisfied  
    # mode = -1 : All cases where x = 1 - μ and dx < 0 are satisfied
    mode = 0

    if mode == 1
        return x[4] >= 0 ? x[1] - (1 - mu) : 1.0
    elseif mode == 0
        return x[1] - (1 - mu)
    elseif mode == -1
        return x[4] < 0 ? x[1] - (1 - mu) : 1.0
    else
        error("Invalid mode: should be -1, 0, or 1")
    end
end

function affect!(stop_order::Int, save_container::Vector{Vector{Float64}})
    """
    This function describes the process to be executed when the condition defined in the above `condition` function is met. 
    Specifically, it saves the state variables and terminates the integration once the specified 
    number of detections has been reached.
    """
    return function (integrator)
        push!(save_container, copy(integrator.u))
        if length(save_container) == stop_order
            terminate!(integrator)
        end
    end
end
