"""
Ex6_multiple_shooting.jl
"""

using DifferentialEquations
using LinearAlgebra
using GLMakie
using Polynomials
using Printf

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

# Retrieving Parameters for the Earth-moon Circular Restricted Three-Body Problem
mu, a_1, a_s, w_1 = fun_cr3bp_parameter(1)
# Calculation of Lagrange points
L1, L2, L3, L4, L5 = fun_libration_points(mu)

"""differential correction"""
iteration_DC_max = 100
threshold = 1e-12

# L1 Lyapunov orbit
x0_1 = [0.987933580858119 0 0 0 0.018052674638026 0]'
t0_1 = 3.290337439058733/2

for i = 1:iteration_DC_max
    global x0_1,t0_1,x_corrected
    global x_n_1,t_n_1,C_L1 = fun_differential_correction_cr3bp(x0_1, t0_1 ,mu)
    global prob_L1 = ODEProblem(fun_cr3bp!, x_n_1, (0.0, 2 * t_n_1), mu)
    global sol_L1 = solve(prob_L1, Vern7(), abstol=3e-14, reltol=1e-14)

    t_corrected = sol_L1.t
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
x0_2 = [1.012164676114102 0 0 0 -0.017947125645568 0]'
t0_2 = 3.319976286115385/2

for i in 1:iteration_DC_max
    global x0_2,t0_2,x_corrected
    global x_n_2,t_n_2,C_L2 = fun_differential_correction_cr3bp(x0_2, t0_2 ,mu)
    global prob_L2 = ODEProblem(fun_cr3bp!, x_n_2, (0.0, 2 * t_n_2), mu)
    global sol_L2 = solve(prob_L2, Vern7(), abstol=3e-14, reltol=1e-14)

    t_corrected = sol_L2.t
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

"""initial trajectory"""
x1 = [0.985000000000000 0.001092112291766 0 0.022 0.005 0]'
x2 = [1-mu -4e-3 0 0.03 0 0]'
x3 = [1.015000000000000 0.001092400176687 0 0.020764610600732 -0.004967672676993 0]'
delta_t1 = 0.75
delta_t2 = 0.75

tspan1 = (0,delta_t1)
prob1 = ODEProblem(fun_cr3bp!, x1, tspan1, mu)
sol1 = solve(prob1, Vern7(), abstol=3e-14, reltol=1e-14)

tspan2 = (0,delta_t2)
prob2 = ODEProblem(fun_cr3bp!, x2, tspan2, mu)
sol2 = solve(prob2, Vern7(), abstol=3e-14, reltol=1e-14)

"""multiple shooting"""
iteration_MS_max = 100
X0 = [x1; x2; x3; delta_t1; delta_t2]
n = 3

for i in 1:iteration_MS_max
    global X0, n
    global X_n = fun_multiple_shooting_cr3bp(X0, n, mu)
    global Xf = zeros(6*(n-1),1)
    global Phi = zeros(6,6,n-1)

    for i = 1:(n-1)
        local Y0 = vcat(X_n[(6*i-5):6*i], reshape(Matrix(I, 6, 6), :))
        local tspan = (0,X_n[6*n+i])
        local prob = ODEProblem(fun_stm_cr3bp!, Y0, tspan, mu)
        global sol = solve(prob, Vern7(), abstol=3e-14, reltol=1e-14)
    
        Xf[(6*i-5):6*i] = sol[1:6,end]'
        Phi[:, :, i] = reshape(sol[7:end,end], 6, 6)
    end

    F = zeros(6*(n-1),1)
    for i = 1:(n-1)
        global F[(6*i-5):6*i] = Xf[(6*i-5):6*i] - X0[(6*i+1):6*(i+1)]
    end

    X_error = norm(F)
    println(@sprintf("Iteration %d: X_error = %.4e", i, X_error))
    
    # The error has converged.
    if X_error < threshold
        break
    end

    # Stop the calculation if it diverges. 
    if X_error > 1e3
        println("calculation diverged in multiple shooting")
        return
    end

    # Within the maximum number of iterations, the error does not fall below the convergence threshold.
    if i == iteration_MS_max
        disp("do not finish in multiple shooting");
        return
    end
    
    # Update X0 
    X0 = X_n
end

fig1 = Figure(size = (800, 600))
ax = Axis(fig1[1,1], xlabel="x [-]", ylabel="y [-]",
        xlabelfont=font, ylabelfont=font,
        xticklabelfont=font, yticklabelfont=font
    )

# Lyapunov orbits around L1 and L2
GLMakie.lines!(ax, sol_L1[1, :], sol_L1[2, :], color=:black, linestyle=:dot)
GLMakie.lines!(ax, sol_L2[1, :], sol_L2[2, :], color=:black, linestyle=:dot)

# Plot Lagrange points
GLMakie.scatter!(ax, [L1[1]], [L1[2]], [L1[3]], color=:black, marker=:star6, markersize=15)
GLMakie.scatter!(ax, [L2[1]], [L2[2]], [L2[3]], color=:black, marker=:star6, markersize=15)

# initial guess1
initial_guess1 = lines!(ax, sol1[1, :], sol1[2, :], color=:blue)
initial_guess2 = lines!(ax, sol2[1, :], sol2[2, :], color=:blue)

Legend(fig1[1, 2], [[initial_guess1]], ["Initial guess"], labelfont = font)

f1_name = "Ex6_initial_guess_mu=$(mu)"
f1_name = replace(f1_name, "." => ",") 
save_path = joinpath(save_fig_dir, f1_name * ".png")
save(save_path, fig1)

fig2 = Figure(size = (800, 600))
ax = Axis(fig2[1,1], xlabel="x[km]", ylabel="y[km]",
        xlabelfont=font, ylabelfont=font,
        xticklabelfont=font, yticklabelfont=font
    )

# Lyapunov Orbits around L1 and L2
GLMakie.lines!(ax, sol_L1[1, :], sol_L1[2, :], color=:black, linestyle=:dot)
GLMakie.lines!(ax, sol_L2[1, :], sol_L2[2, :], color=:black, linestyle=:dot)

# Plot Lagrange points
GLMakie.scatter!(ax, [L1[1]], [L1[2]], [L1[3]], color=:black, marker=:star6, markersize=15)
GLMakie.scatter!(ax, [L2[1]], [L2[2]], [L2[3]], color=:black, marker=:star6, markersize=15)

multiple_shooting_plots = []

for i = 1:(n-1)
    local tspan = (0,X_n[6*n+i])
    local prob = ODEProblem(fun_cr3bp!, X_n[(6*i-5):6*i], tspan, mu)
    local sol = solve(prob, Vern7(), abstol=3e-14, reltol=1e-14)

    push!(multiple_shooting_plots, lines!(ax, sol[1,:], sol[2,:], color = "#D95319"))

    GLMakie.scatter!(ax, [sol[1,1]], [sol[2,1]], color = "#D95319", markersize = 10)
    GLMakie.scatter!(ax, [sol[1,end]], [sol[2,end]], color = "#D95319", markersize = 10)
end

# Initial guess 
initial_guess1 = lines!(ax, sol1[1, :], sol1[2, :], color=:blue)
initial_guess2 = lines!(ax, sol2[1, :], sol2[2, :], color=:blue)

Legend(fig2[1, 2], [[initial_guess1], [multiple_shooting_plots[1]]], 
    ["Initial guess","Multiple shooting"], labelfont = font)

f2_name = "Ex6_multiple_shooting_mu=$(mu)"
f2_name = replace(f2_name, "." => ",") 
save_path = joinpath(save_fig_dir, f2_name * ".png")
save(save_path, fig2)

display(fig2)