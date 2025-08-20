"""
Ex4_differential_correction.jl
"""

using DifferentialEquations
using LinearAlgebra
using GLMakie
using Polynomials
using Printf

# Function loading
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
# Save Figure directory
save_fig_dir = "Figure/"
# Save dat directory
save_dat_dir = "Output_data/"

# Retrieving Parameters for the Sun-Earth Circular Restricted Three-Body Problem
mu, a_1, a_s, w_1 = fun_cr3bp_parameter(1)
# Calculation of Lagrange points
L1, L2, L3, L4, L5 = fun_libration_points(mu)

# Initial conditions for Lyapunov orbit
x0 = [0.9889 0 0 0 0.008375 0]'
t0 = 2

# Initial conditions for halo orbit
#x0 = [0.9919, 0, 0.002223, 0, -0.01014, 0]
#t0 = pi / 2

if x0[3] == 0
	Lyapunov = true;
else
	Lyapunov = false;
end

"""differential_correction"""
parameter = (mu)
# Define ODE Problem
tspan = (0.0, 2 * t0)

# Initial guess
prob = ODEProblem(fun_cr3bp!, x0, tspan, parameter)
sol_init = solve(prob, Vern7(), abstol=3e-14, reltol=1e-14)

# Maximum number of iterations for differential correction
iteration_max = 100
# Convergence threshold
threshold = 1e-11

# Correction
for i in 1:iteration_max
    global x0,t0,x_corrected
    global x_n, t_n, C = fun_differential_correction_cr3bp(x0,t0,mu)
    global prob = ODEProblem(fun_cr3bp!, x0, (0.0, 2 * t_n), mu)
    global sol = solve(prob, Vern7(), abstol=3e-14, reltol=1e-14)

    t_corrected = sol.t
    x_corrected = sol.u

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
    if i == iteration_max
        println("do not finish");
        return
    end
    
    # Update x0 and t0 
    x0 = x_n
    t0 = t_n
end

# Figure
fig1 = Figure(size=(800, 600))

if Lyapunov 
    ax = Axis3(fig1[1,1], xlabel="x[-]", ylabel="y[-]", zlabel=L"z[-]", azimuth = -pi/2, elevation = pi/2, 
        xlabelfont = font, ylabelfont = font, zlabelfont = font, xticklabelfont = font, yticklabelfont = font, zticklabelfont = font
    )
    ax.zlabelvisible = false
    ax.zticklabelsvisible = false
  else
    ax = Axis3(fig1[1,1], xlabel=L"x[-]", ylabel=L"y[-]", zlabel=L"z[-]", azimuth = pi/9, elevation = pi/6,
        xlabelfont = font, ylabelfont = font, zlabelfont = font, xticklabelfont = font, yticklabelfont = font, zticklabelfont = font
    )
end

# Initial orbit (red line)
lines!(ax, sol_init[1,:], sol_init[2,:], sol_init[3,:], color=:red)

# Plot initial point of the initial orbit (red dot)
iniial = GLMakie.scatter!(ax, [sol_init[1,1]], [sol_init[2,1]], [sol_init[3,1]], 
         markersize=10, color=:red, marker=:circle)

# Periodic orbit (blue line)
lines!(ax, sol[1,:], sol[2,:], sol[3,:], color=:blue)

# Plot initial point of the periodic orbit (blue dot)
corrected = GLMakie.scatter!(ax, [sol[1,1]], [sol[2,1]], [sol[3,1]], 
         markersize=10, color=:blue, marker=:circle)

# Plot Lagrange points
L1_point = GLMakie.scatter!(ax, [L1[1]], [L1[2]], [L1[3]], color=:black, marker=:star6, markersize=15)

Legend(fig1[1, 2], [iniial, corrected, L1_point], ["iniial", "corrected", "L1"], labelfont = font)

f1_name = "Ex4_differential_correction_mu=$(mu)_x0=$(x0[1])_z0=$(x0[3])_ydot0=$(x0[5])_t0=$(t0)"
f1_name = replace(f1_name, "." => ",")  
save_path = joinpath(save_fig_dir, f1_name * ".png")  
save(save_path, fig1)  

fun_save_data(sol.t, sol.u, save_dat_dir, f1_name)

display(fig1)  

"""monodromy matrix"""
X0 = vcat(x_n, reshape(Matrix(I, 6, 6), :))
tspan = (0.0, 2 * t_n)

prob = ODEProblem(fun_stm_cr3bp!, X0, tspan, parameter)
sol = solve(prob, Vern7(), abstol=3e-14, reltol=1e-14)

monodromy = reshape(sol.u[end][7:end], 6, 6)
D, V = eigen(monodromy)

# Data for eigenvalue plotting
theta = range(0, 2π, length=1001)
x_circle = cos.(theta)
y_circle = sin.(theta)

# Retrieve the eigenvalues of the monodromy matrix
eigenvalues = D

# Plot of the unit circle and eigenvalues
fig2 = Figure()
ax1 = Axis(fig2[1, 1], xlabel="real part", ylabel="imaginary part", aspect=1,
    xlabelfont = font, ylabelfont = font, xticklabelfont = font, yticklabelfont = font
)
lines!(ax1, x_circle, y_circle, color=:black)

scatter_plot = GLMakie.scatter!(ax1, real.(eigenvalues), imag.(eigenvalues), color=:blue, markersize=10, marker=:circle)
Legend(fig2[1, 2], [scatter_plot], ["eigenvalue"], labelfont = font)

f2_name = "Ex4_monodromy_mu=$(mu)_x0=$(x0[1])_z0=$(x0[3])_ydot0=$(x0[5])_t0=$(t0)"
f2_name = replace(f2_name, "." => ",")
save_path = joinpath(save_fig_dir, f2_name * ".png")
save(save_path, fig2)

display(fig2)

# Plot with restricted coordinate range
fig3 = Figure()
ax2 = Axis(fig3[1, 1], xlabel="real part", ylabel="imaginary part", aspect=1, limits=(-1.5, 1.5, -1.5, 1.5), 
    xlabelfont = font, ylabelfont = font, xticklabelfont = font, yticklabelfont = font
)
lines!(ax2, x_circle, y_circle, color=:black)
GLMakie.scatter!(ax2, real.(eigenvalues), imag.(eigenvalues), color=:blue, markersize=10, marker=:circle)
Legend(fig3[1, 2], [scatter_plot], ["eigenvalue"], labelfont = font)

f3_name = "Ex4_monodromy_origin_mu=$(mu)_x0=$(x0[1])_z0=$(x0[3])_ydot0=$(x0[5])_t0=$(t0)"
f3_name = replace(f3_name, "." => ",")
save_path = joinpath(save_fig_dir, f3_name * ".png")
save(save_path, fig3)

display(fig3)