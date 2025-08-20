"""
Ex8_poincare_map.jl
"""

using DifferentialEquations
using LinearAlgebra
using GLMakie
using Polynomials
using Printf
using CairoMakie
using Colors

current_dir = pwd()

# Function loading
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

# Retrieving Parameters for the Sun-Jupiter Circular Restricted Three-Body Problem
mu, a_1, a_s, w_1 = fun_cr3bp_parameter(6)
# Calculation of Lagrange points
L1, L2, L3, L4, L5 = fun_libration_points(mu)

"""differential correction"""
iteration_DC_max = 100;
threshold = 1e-13;

# L1 Lyapunov orbit
x0_1 = [0.926672517344211 0 0  0 0.045684284381193 0]'
t0_1 = 2.915640890603734/2

for i = 1:iteration_DC_max
    global x0_1,t0_1,x_corrected
    global x_n_1,t_n_1,C_L1 = fun_differential_correction_cr3bp(x0_1, t0_1 ,mu)
    global prob_L1 = ODEProblem(fun_cr3bp!, x_n_1, (0.0, 2 * t_n_1), mu)
    global sol_L1 = solve(prob_L1, Vern7(), abstol=3e-14, reltol=1e-14)

    x_corrected = sol_L1.u

    # Calculate the difference between the initial state and the state after one orbital period.
    x_error = norm(x_corrected[end] .- x_corrected[1], 2)
    println(@sprintf("Iteration %d: x_error = %.4e", i, x_error))
    
    # The error has converged.
    if x_error < threshold
        break
    end

    # Stop the calculation if it diverges. 
    if x_error > 1e3
        println("Calculation diverged")
        return
    end

    # Within the maximum number of iterations, the error does not fall below the convergence threshold.
    if i == iteration_DC_max
        println("do not finish");
        return
    end
    
    # Update x0 and t0 
    x0_1 = x_n_1
    t0_1 = t_n_1
end

x_L1 = x_corrected
t_L1 = t_n_1

# L2 Lyapunov orbit
x0_2 = [1.072915577273646 0 0 0 -0.025945940908770 0]'
t0_2 = 3.184772418214781/2

for i in 1:iteration_DC_max
    global x0_2,t0_2,x_corrected
    global x_n_2,t_n_2,C_L2 = fun_differential_correction_cr3bp(x0_2, t0_2 ,mu)
    global prob_L2 = ODEProblem(fun_cr3bp!, x_n_2, (0.0, 2 * t_n_2), mu)
    global sol_L2 = solve(prob_L2, Vern7(), abstol=3e-14, reltol=1e-14)

    x_corrected = sol_L2.u

    # Calculate the difference between the initial state and the state after one orbital period.
    x_error = norm(x_corrected[end] .- x_corrected[1], 2)
    println(@sprintf("Iteration %d: x_error = %.4e", i, x_error))
    
    # The error has converged.
    if x_error < threshold
        break
    end

    # Stop the calculation if it diverges. 
    if x_error > 1e3
        println("Calculation diverged")
        return
    end

    # Within the maximum number of iterations, the error does not fall below the convergence threshold.
    if i == iteration_DC_max
        println("do not finish");
        return
    end
    
    # Update x0 and t0 
    x0_2 = x_n_2
    t0_2 = t_n_2
end

x_L2 = x_corrected
t_L2 = t_n_2

"""zero_velocity_curve"""
x_range = 0.9:1e-4:1.1
y_range = -0.1:1e-4:0.1
x = repeat(collect(x_range), 1, length(y_range)) 
y = repeat(collect(y_range)', length(x_range), 1)  
z = zeros(size(x))  

# Potential function calculation
r1 = sqrt.((x .+ mu) .^ 2 .+ y .^ 2 .+ z .^ 2)
r2 = sqrt.((x .- (1 - mu)) .^ 2 .+ y .^ 2 .+ z .^ 2)
U = 0.5 .* (x .^ 2 .+ y .^ 2) .+ (1 - mu) ./ r1 .+ mu ./ r2
C = 2 .* U

tspan_u = (0,9)
tspan_s = (10,0)
xpert = 1e-10
N = 40

# Initialize storage arrays:
yu_right_f1 = zeros(N, 6)
yu_right_f2 = fill(NaN, N, 6)
ys_left_f1  = zeros(N, 6)
ys_left_f2  = fill(NaN, N, 6)

# L1 unstable manifold generation
_, _, _, XU_right, _ = fun_manifold_cr3bp(mu, x_L1[1], 2 * t_L1, N, xpert)
# L2 stable manifold generation
XS_left, _, _, _, _ = fun_manifold_cr3bp(mu, x_L2[1], 2 * t_L2, N, xpert)

# Create a new figure and add a filled contour plot.
fig = Figure(size=(600, 600))

# Create an Axis for the Plot
ax = Axis(fig[1, 1]; aspect = 1, xlabel = "x[-]", ylabel = "y [-]",
        xlabelfont=font, ylabelfont=font,
        xticklabelfont=font, yticklabelfont=font
    )

# Create a filled contour at the level C_L1
GLMakie.contourf!(ax, x_range, y_range, C; levels=[minimum(C), C_L1], colormap=[:gray])
GLMakie.contour!(ax, x_range, y_range, C; levels = [C_L1], colormap = :grays)

# Right branch integration
for i in 1:N
    # Create an empty vector to store event states for the current integration
    local xe = Vector{Vector{Float64}}()

    # define a callback
    local cb_local = ContinuousCallback(condition, affect!(3, xe); save_positions=(false, false))
    
    # Define the ODE problem.
    local parameter = (mu)
    local prob = ODEProblem(fun_cr3bp!, XU_right[:, i], tspan_u, parameter)
    
    # Solve the ODE with high tolerances and the callback.
    local sol = solve(prob, Vern7(); reltol=3e-14, abstol=1e-14, callback=cb_local)
    # Reverse the trajectory (flip vertically) as in MATLAB's flipud.
    
    # Plot the reversed trajectory in color "#77AC30".
    GLMakie.lines!(sol[1,:], sol[2,:], color="red")
    
    # Save the event states.
    if length(xe) >= 1
        yu_right_f1[i, :] = xe[1]
    end
    if length(xe) >= 2
        yu_right_f2[i, :] = xe[2]
    end 
end

# --- Left branch integration ---
for i in 1:N
    local xe = Vector{Vector{Float64}}()

    local cb_local = ContinuousCallback(condition, affect!(3, xe); save_positions=(false, false))
    
    # Define the ODE problem.
    local parameter = (mu)
    local prob = ODEProblem(fun_cr3bp!, XS_left[:, i], tspan_s, parameter)
    
    # Solve the ODE with high tolerances and the callback.
    local sol = solve(prob, Vern7(); reltol=3e-14, abstol=1e-14, callback=cb_local)
    # Reverse the trajectory (flip vertically) as in MATLAB's flipud.
    
    # Plot the reversed trajectory in color "#77AC30".
    GLMakie.lines!(sol[1,:], sol[2,:], color="green")
    
    # Save the event states.
    if length(xe) >= 1
        ys_left_f1[i, :] = xe[1]
    end
    if length(xe) >= 2
        ys_left_f2[i, :] = xe[2]
    end
end

xlims!(ax, 0.9, 1.1)
ylims!(ax, -0.08, 0.08)

ax.xgridvisible[] = true
ax.ygridvisible[] = true

# Plot unstable manifold (x_L1) in black and assign a label.
p1 = GLMakie.lines!(ax, [v[1] for v in x_L1], [v[2] for v in x_L1], color = :black, label = "Unstable manifold")
# Plot stable manifold (x_L2) in black and assign a label.
p2 = GLMakie.lines!(ax, [v[1] for v in x_L2], [v[2] for v in x_L2], color = :black, label = "Stable manifold")

# Plot Lagrange points as star markers
GLMakie.scatter!(ax, [L1[1]], [L1[2]], markersize = 15, marker = :star5, color = :black)
GLMakie.scatter!(ax, [L2[1]], [L2[2]], markersize = 15, marker = :star5, color = :black)

p1 = GLMakie.lines!(ax, [v[1] for v in x_L1], [v[2] for v in x_L1], color = :black, label = "Unstable manifold")
# Plot stable manifold (x_L2) in black and assign a label.
p2 = GLMakie.lines!(ax, [v[1] for v in x_L2], [v[2] for v in x_L2], color = :black, label = "Stable manifold")

# Plot Lagrange points as star markers
GLMakie.scatter!(ax, [L1[1]], [L1[2]], markersize = 15, marker = :star5, color = :black)
GLMakie.scatter!(ax, [L2[1]], [L2[2]], markersize = 15, marker = :star5, color = :black)

# Add a vertical line at x = 1 - mu with a dashed line style.
GLMakie.lines!(ax, [1 - mu, 1 - mu], [-0.08, 0.08], linestyle = :dash, color = :black, linewidth = 1.5)

# Create a name for saving the figure
f1_name = "Ex8_poincare_map_traj_mu=$(mu)_xpert=$(xpert)"
f1_name = replace(f1_name, "." => ",")
# Display figure 1
display(fig)

# Poincare Map Plot
f2 = Figure(size = (800, 600))
ax2 = Axis(f2[1, 1], xlabel ="y[-]", ylabel ="y_dot[-]",
        xlabelfont=font, ylabelfont=font,
        xticklabelfont=font, yticklabelfont=font
    )

# Define a color for stable manifold insertion (hex color "#77AC30")
color_stable = parse(Colorant, "#77AC30")

# Plot the four insertion curves:
p2_1 = lines!(ax2, yu_right_f1[:, 2], yu_right_f1[:, 5], color = :red, label = "1st insertion of unstable manifold")
p2_2 = lines!(ax2, ys_left_f1[:, 2], ys_left_f1[:, 5], color = color_stable, label = "1st insertion of stable manifold")
p2_3 = lines!(ax2, yu_right_f2[:, 2], yu_right_f2[:, 5], color = :magenta, label = "2nd insertion of unstable manifold")
p2_4 = lines!(ax2, ys_left_f2[:, 2], ys_left_f2[:, 5], color = :green, label = "2nd insertion of stable manifold")

# Set axis square by matching data aspect ratio 
leg2 = Legend(f2, [p2_1, p2_2, p2_3, p2_4],
    ["1st insertion of unstable manifold", "1st insertion of stable manifold",
     "2nd insertion of unstable manifold", "2nd insertion of stable manifold"],
    tellwidth = false, markersize = 10, labelfont = font
)
# You can adjust the legend position by modifying layout; here we add it to the right of the axis.
f2[1, 2] = leg2

# Create a name for saving the figure.
f2_name = "Ex8_poincare_map_mu=$(mu)_xpert=$(xpert)"
f2_name = replace(f2_name, "." => ",")

# Display figure 2
display(f2)