"""
Ex5_manifold.jl
"""

using DifferentialEquations
using LinearAlgebra
using GLMakie
using Polynomials
using Printf

GLMakie.activate!()
# Function loading
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

# Retrieving Parameters for the Sun-Earth Circular Restricted Three-Body Problem
mu, a_1, a_s, w_1 = fun_cr3bp_parameter(1)
# Calculation of Lagrange points
L1, L2, L3, L4, L5 = fun_libration_points(mu)

# Initial conditions for Lyapunov orbit
x0 = [0.9889 0 0 0 0.008375 0]'
#x0 = [0.9919 0 0.002223 0 -0.01014 0]'
t0 = 2

if x0[3] == 0
	Lyapunov = true;
else
	Lyapunov = false;
end

"""differential_correction"""
# Parameter setting
parameter = (mu)
# Define ODE Problem
tspan = (0.0, 2 * t0)

# Initial guess
prob = ODEProblem(fun_cr3bp!, x0, tspan, parameter)
sol_init = solve(prob, Vern7(), abstol=3e-14, reltol=1e-14)

# Maximum number of iterations for differential correction
iteration_max = 100
# Convergence threshold
threshold = 1e-12

# Correction
for i in 1:iteration_max
    global x0,t0,x_corrected
    global x_n, t_n, C_xn = fun_differential_correction_cr3bp(x0,t0,mu)
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

"""Zero_velocity_curve"""
# Mesh generation
x_range =  0.9:1e-4:1.1
y_range = -0.1:1e-4:0.1
x = repeat(collect(x_range), 1, length(y_range)) 
y = repeat(collect(y_range)', length(x_range), 1)  
z = zeros(size(x))  

# Potential function calculation
r1 = sqrt.((x .+ mu) .^ 2 .+ y .^ 2 .+ z .^ 2)
r2 = sqrt.((x .- (1 - mu)) .^ 2 .+ y .^ 2 .+ z .^ 2)
U = 0.5 .* (x .^ 2 .+ y .^ 2) .+ (1 - mu) ./ r1 .+ mu ./ r2
C = 2 .* U

"""stable and unstable manifolds"""
xpert = 1e-6
N = 40

XS_left, XS_right, XU_left, XU_right, Y = fun_manifold_cr3bp(mu, x_n, 2*t_n, N, xpert)

tf = 4.8

tspan_s = (tf, 0)
tspan_u = (0, tf)

fig = Figure(size=(800, 600))
if Lyapunov
    ax = Axis(fig[1, 1], xlabel="x[-]", ylabel="y[-]",
        xlabelfont=font, ylabelfont=font,
        xticklabelfont=font, yticklabelfont=font
    )
else
    ax = Axis3(fig[1, 1], xlabel="x[-]", ylabel="y[-]", zlabel=L"z[-]",
        azimuth=pi/9, elevation=pi/6,
        xlabelfont=font, ylabelfont=font, zlabelfont=font,
        xticklabelfont=font, yticklabelfont=font, zticklabelfont=font
    )
end

# Stable manifolds
f1_p1 = []
f1_p2 = []
for i in 1:N
    prob_ys_left = ODEProblem(fun_cr3bp!, XS_left[:, i], (tspan_s[1], tspan_s[end]), parameter)
    sol_ys_left = solve(prob_ys_left, Vern7(), abstol=1e-14, reltol=1e-14)
    if Lyapunov
        push!(f1_p1, lines!(ax, sol_ys_left[1, :], sol_ys_left[2, :], color=:green))
    else
        push!(f1_p1, lines!(ax, sol_ys_left[1, :], sol_ys_left[2, :], sol_ys_left[3, :], color=:green))
    end
end

for i in 1:N
    prob_ys_right = ODEProblem(fun_cr3bp!, XS_right[:, i], (tspan_s[1], tspan_s[end]), parameter)
    sol_ys_right = solve(prob_ys_right, Vern7(), abstol=1e-14, reltol=1e-14)
    if Lyapunov
        push!(f1_p2, lines!(ax, sol_ys_right[1,:], sol_ys_right[2,:], color=:blue))
    else
        push!(f1_p2, lines!(ax, sol_ys_right[1,:], sol_ys_right[2,:], sol_ys_right[3,:], color=:blue))
    end
end

# Unstable manifolds
f1_p3 = []
f1_p4 = []
for i in 1:N
    prob_yu_left = ODEProblem(fun_cr3bp!, XU_left[:, i], (tspan_u[1], tspan_u[end]), parameter)
    sol_yu_left = solve(prob_yu_left, Vern7(), abstol=1e-14, reltol=1e-14)
    if Lyapunov
        push!(f1_p3, lines!(ax, sol_yu_left[1,:], sol_yu_left[2,:], color=:magenta))
    else
        push!(f1_p3, lines!(ax, sol_yu_left[1,:], sol_yu_left[2,:], sol_yu_left[3,:], color=:magenta))
    end
end

for i in 1:N
    prob_yu_right = ODEProblem(fun_cr3bp!, XU_right[:, i], (tspan_u[1], tspan_u[end]), parameter)
    sol_yu_right = solve(prob_yu_right, Vern7(), abstol=1e-14, reltol=1e-14)
    if Lyapunov
        push!(f1_p4, lines!(ax, sol_yu_right[1,:], sol_yu_right[2,:], color=:red))
    else
        push!(f1_p4, lines!(ax, sol_yu_right[1,:], sol_yu_right[2,:], sol_yu_right[3,:], color=:red))
    end
end

# Legend
Legend(fig[1, 2], [f1_p1[1], f1_p2[1], f1_p3[1], f1_p4[1]], 
["left half of stable manifold", "right half of stable manifold", "left half of unstable manifold", "right half of unstable manifold"]
, labelfont = font)

# Plot corrected trajectory and L1 point
if Lyapunov
    lines!(ax, sol[1, :], sol[2, :], color=:black)
    GLMakie.scatter!(ax, [L1[1]], [L1[2]], markersize=10, color=:black)
    xlims!(ax, 0.96, 1.02)
    ylims!(ax, -0.02, 0.02)
else
    lines!(ax, sol[1, :], sol[2, :], sol[3, :], color=:black)
    GLMakie.scatter!(ax, [L1[1]], [L1[2]], [L1[3]], markersize=10, color=:black)
    xlims!(ax, 0.96, 1.02)
    ylims!(ax, -0.02, 0.02)
    zlims!(ax, -0.02, 0.02)
end

# Save figure
f1_name = "Ex5_manifold_mu=$(mu)_xn=$(x_n[1])_zn=$(x_n[3])_ydotn=$(x_n[5])_tn=$(t_n)_xpert=$(xpert)_tf=$(tf)"
f1_name = replace(f1_name, "." => ",")

if Lyapunov 
    # Zero velocity curve
    contourf!(ax, x_range, y_range, C; levels=[minimum(C), C_xn], colormap=[:gray])
    contour!(ax, x_range, y_range, C; levels=[C_xn], linewidth=1.5, color=:black)
    save_path = joinpath(save_fig_dir, f1_name * ".png")
    save(save_path, fig)
else
    ax.azimuth = pi/2
    ax.elevation = pi
    ax.ylabel = ""
    ax.yticklabelsvisible=false
    save_path = joinpath(save_fig_dir, f1_name * "_xz.png")
    save(save_path, fig)
    ax.azimuth = 0.0
    ax.elevation = 0.0
    ax.ylabel = "y"
    ax.yticklabelsvisible=true
    ax.xlabel = ""
    ax.xticklabelsvisible=true
    save_path = joinpath(save_fig_dir, f1_name * "_yz.png")
    save(save_path, fig)
    # Zero velocity curve
    contourf!(ax, x_range, y_range, C; levels=[minimum(C), C_xn], colormap=[:gray])
    contour!(ax, x_range, y_range, C; levels=[C_xn], linewidth=1.5, color=:black)
    ax.azimuth = -pi/2
    ax.elevation = pi/2
    ax.xlabel = "x"
    ax.xticklabelsvisible=true
    ax.zlabel = ""
    ax.zticklabelsvisible=false
    save_path = joinpath(save_fig_dir, f1_name * "_xy.png")
    save(save_path, fig)
end

display(fig)  