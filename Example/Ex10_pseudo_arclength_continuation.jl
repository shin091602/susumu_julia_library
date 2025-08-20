"""
Ex10_pseudo_arclength_continuation.jl
"""

using GLMakie
using CairoMakie
using DifferentialEquations
using Polynomials
using ColorSchemes
using Printf
using LinearAlgebra
using Interpolations 

GLMakie.activate!()
#GLMakie.inline!(false)

# Load functions
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

# Retrieve parameters for the Earth-Moon Circular Restricted Three-Body Problem
mu, a_1, a_s, w_1 = fun_cr3bp_parameter(1)
dim_L = a_1
# Compute Lagrange points
L1, L2, L3, L4, L5 = fun_libration_points(mu)

# Iteration and plot settings
# Max number of iterations for differential correction
iteration_max = 100          
# Convergence threshold for correction
threshold     = 1e-10
# Step scale for pseudo-arclength continuation
scale         = 1e-4 
# Max number of total steps
count_max     = 400
# Plot update interval
plot_interval = 100

# Jacobi constant range (for colormap)
Jacobi_min = 3.000212571920089
Jacobi_max = 3.000818466710653

# Plotting range (axis limits)
xlim_min = 1.485e8
xlim_max = 1.525e8
ylim_min = -2.5e6
ylim_max =  2.5e6
zlim_min = -2.3e6
zlim_max =  1.0e6

# Setup colormap 
# Create 256-color "jet" colormap and assign it to Jacobi constant range
color_array = collect(cgrad(:jet, 256, categorical=false))
Jacobi_lim = range(Jacobi_min, Jacobi_max, length=length(color_array))
Interp_c = interpolate((Jacobi_lim,), color_array, Gridded(Linear()))

# Arrays for storing results
num_points = div(count_max, plot_interval) + 1
Jacobi_arr   = zeros(Float64, num_points)
x0_corrected = zeros(Float64, num_points, 6)
t0_corrected = zeros(Float64, num_points)

# Initial condition
x01_ast = [1.00837693674644, 0.0, 0.000200000000000000, 0.0, 0.00976542047036519, 0.0]
t01_ast = 1.55123023270505

"""pseudo arclength continuation"""
Jacobi_arr[1] = fun_Jacobi_const(x01_ast, mu)
prob = ODEProblem(fun_cr3bp!, x01_ast, (0.0, 2*t01_ast), mu)
sol = solve(prob, Vern7(), abstol=1e-14, reltol=1e-14)
x0_corrected[1, :] = x01_ast'
t0_corrected[1]   = t01_ast

# Create figure and 3D axis 
fig = Figure(size = (1000, 600))
ax = Axis3(fig[1,1], aspect = (1,1,1), xlabel="x[km]", ylabel="y[km]", zlabel="z[km]",
        xlabelfont=font, ylabelfont=font, xticklabelfont=font, yticklabelfont=font)

ax.xticks = 1.48e8:0.01e8:1.52e8

rgb = Interp_c[Jacobi_arr[1]]
GLMakie.lines!(ax, sol[1,:]*dim_L, sol[2,:]*dim_L, sol[3,:]*dim_L, color = rgb)

for i in 1:count_max
    global x01_ast, t01_ast
    local delta, ~ = fun_null_cr3bp(x01_ast, t01_ast, mu)
    delta = -delta
    local x02 = copy(x01_ast) 
    x02[1] = x01_ast[1] + scale * delta[1]
    x02[3] = x01_ast[3] + scale * delta[2]
    x02[5] = x01_ast[5] + scale * delta[3]
    t02    = t01_ast    + scale * delta[4]

    # Differential correction loop (inner iteration until convergence)
    for j in 1:iteration_max
        # Get corrected initial state, period, and Jacobi constant C from correction function
        global x02_ast, t02_ast, C, G = fun_differential_correction_cr3bp_PAC(x01_ast, t01_ast, x02, t02, scale, delta, mu)
        global prob = ODEProblem(fun_cr3bp!, x02_ast, (0.0, 2*t02_ast), mu)
        global sol = solve(prob, Vern7(), abstol=1e-14, reltol=1e-14)

        G_error = norm(G)

        # Exit inner loop if converged
        if G_error < threshold
            break
        end
        
        # Abort if divergence detected
        if G_error > 1e3
            println("Calculation diverged")
            break
        end

        if j == iteration_max
            println("do not finish");
            break
        end
        
        # Update initial condition with current correction
        x02 = x02_ast
		t02 = t02_ast
    end

    # Plot update: draw orbit every plot_interval
    if i >= plot_interval && i % plot_interval == 0
        idx = div(i, plot_interval)
        Jacobi_arr[idx] = C
        x0_corrected[idx, :] = sol.u[1]
        t0_corrected[idx] = sol.t[end]
        
        # Get color from interpolator corresponding to Jacobi constant C
        local rgb = Interp_c[C]
        
        GLMakie.lines!(ax, sol[1,:]*dim_L, sol[2,:]*dim_L, sol[3,:]*dim_L, color = rgb)
    end

    x01_ast = x02_ast
    t01_ast = t02_ast

    println(@sprintf("Iteration %d", i))
end

# Supplementary plot: L₂ point and Earth
# Plot L₂ point (star marker)
f1_p1 = scatter!(ax, [L2[1]*dim_L], [L2[2]*dim_L], [L2[3]*dim_L],
    markersize = 14,
    marker = :star5,
    color = :black
)
# Plot Earth (circle marker)
f1_p2 = scatter!(ax, [(1-mu)*dim_L], [0], [0],
    markersize = 14,
    marker = :circle,
    color = :black
)

# Set view angle, axis limits, and grid
xlims!(ax, (xlim_min, xlim_max))
ylims!(ax, (ylim_min, ylim_max))
zlims!(ax, (zlim_min, zlim_max))

cbar = Colorbar(fig[1, 2],colormap = :jet,limits = (Jacobi_min, Jacobi_max),
    label = "Jacobi constant [-]", labelsize = 15,labelfont = font, ticklabelfont = font
    )
cbar.ticks = LinRange(Jacobi_min, Jacobi_max, 6)

leg = axislegend(ax, [f1_p1, f1_p2], ["L_2", "Earth"], position = :lt, labelfont = font)

display(fig)

f1_name = "Ex10_mu=$(round(mu, sigdigits=3))" *
          "_x=$(round(x0_corrected[1,1], digits=4))" *
          "_z=$(round(x0_corrected[1,3], digits=4))" *
          "_v=$(round(x0_corrected[1,5], digits=4))" *
          "_t=$(round(t0_corrected[1], digits=4))" *
          "_s=$(scale)" *
          "_n=$(count_max)" *
          "_pi=$(plot_interval)"


f1_name = replace(f1_name, "." => ",") 
save_path = joinpath(save_fig_dir, f1_name*".png")
save(save_path, fig)
