"""
Ex9_natural_parameter_continuation.jl
"""

using DifferentialEquations
using LinearAlgebra
using GLMakie
using Polynomials
using Printf
using Interpolations

current_dir = pwd()

functions_dir = normpath(joinpath(current_dir, "..", "Functions"))

if isdir(functions_dir)
    for file in readdir(functions_dir)
        if endswith(file, ".jl")
            full_path = joinpath(functions_dir, file)
            if isfile(full_path)
                include(full_path)
            else
                println("Warning: File not found - ", full_path)
            end
        end
    end
else
    println("Error: Functions directory not found at ", functions_dir)
end

# Start timer
t0 = time()  
# font setting
font = "Times New Roman"
# Save directory
save_fig_dir = "Figure/"

# Retrieving parameters for the Earth-Moon Circular Restricted Three-Body Problem
mu, a_1, a_s, w_1 = fun_cr3bp_parameter(1)
dim_L = a_1
# Compute Lagrange points
L1, L2, L3, L4, L5 = fun_libration_points(mu)

iteration_max = 100
threshold = 1e-10
delta = 2e-6
count = 0
# Max number of total steps
count_max = 4000
# Plot update interval
plot_interval = 200

# Jacobi constant range (for colormap)
Jacobi_min = 3.000040465479221
Jacobi_max = 3.000881506800496

# Maximum number of iterations for differential correction
iteration_max = 100
# Convergence threshold
threshold = 1e-10

color_grad = cgrad(:jet)

Jacobi = zeros(Float64, count_max ÷ plot_interval)
x0_corrected = zeros(Float64, count_max ÷ plot_interval, 6)
t0_corrected = zeros(Float64, count_max ÷ plot_interval)

x0 = [1.01003168694194, 0, 0, 0, 8.61411407047802e-6, 0]
t0 = 1.47271346920646

fig = Figure(size = (800, 600))
ax = Axis(fig[1, 1], xlabel="x[km]", ylabel="y[km]", aspect = DataAspect(), xlabelfont=font, ylabelfont=font,
        xticklabelfont=font, yticklabelfont=font)

while true
    global count
    count += 1

    for i in 1:iteration_max
        global x0, t0, x_corrected
        global x_n, t_n, C = fun_differential_correction_cr3bp(x0, t0, mu)
        global prob = ODEProblem(fun_cr3bp!, x0, (0.0, 2 * t_n), mu)
        global sol = solve(prob, Vern7(), abstol=3e-14, reltol=1e-14)

        # Calculate the difference between the initial state and the state after one orbital period
        x_error = norm(sol.u[end] .- sol.u[1], 2)
        println(@sprintf("Iteration %d: x_error = %.4e", count, x_error))
        
        # Exit loop if convergence criterion is met
        if x_error < threshold
            break
        end
        
        # Abort if the solution diverges
        if x_error > 1e3
            println("Calculation diverged")
            return
        end

        if i == iteration_max
            println("do not finish")
            return
        end
        
        # Update x0 and t0
        x0 = x_n
        t0 = t_n
    end

    # Update plot every plot_interval steps
    if count >= plot_interval && count % plot_interval == 0
        idx = div(count, plot_interval)
        Jacobi[idx] = C
        x0_corrected[idx, :] = sol.u[1]
        t0_corrected[idx] = sol.t[end]

        # Normalize C into the range [0, 1]
        norm_val = (C - Jacobi_min) / (Jacobi_max - Jacobi_min)
        local rgb = color_grad[norm_val]  

        lines!(ax, sol[1, :] * dim_L, sol[2, :] * dim_L, color = rgb)
    end
    
    # Decrease the first element of x0 by delta
    x0[1] -= delta

    if count == count_max
        break
    end
end

# Plot L2 point and Earth
scatter!(ax, [L2[1] * dim_L], [L2[2] * dim_L], color = :black, marker = :star6, markersize = 15)
scatter!(ax, [(1 - mu) * dim_L], [0], color = :black, marker = :star6, markersize = 15)

# Add colorbar
cb = Colorbar(fig[1, 2], colormap = color_grad, limits = (Jacobi_min, Jacobi_max),
              label = "Jacobi constant [-]", width = 20, labelfont = font, ticklabelfont = font
    )

#f1_name = "Ex9_natural_parameter_continuation_mu=$(mu)_x0=$(x0[1])_ydot0=$(x0[5])_t0=$(t0)_threshold=$(threshold)_delta=$(delta)_count_max=$(count_max)_plot_interval=$(plot_interval)"
f1_name = "Ex9_natural_parameter_continuation_mu=$(mu)_x0=$(x0[1])_ydot0=$(x0[5])_t0=$(t0)_threshold=$(threshold)"
f1_name = replace(f1_name, "." => ",")
save_path = joinpath(save_fig_dir, f1_name * ".png")
save(save_path, fig)

display(fig)
