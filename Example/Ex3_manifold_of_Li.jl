"""
Ex3_manifold_of_Li.jl
"""

using DifferentialEquations
using LinearAlgebra
using GLMakie
using Polynomials

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

"""Zero_velocity_curve"""
# Retrieving Parameters for the Earth-Moon Circular Restricted Three-Body Problem
mu, a_1, a_s, w_1 = fun_cr3bp_parameter(2)
# Calculation of Lagrange points
L1, L2, L3, L4, L5 = fun_libration_points(mu)

# Mesh generation
x_range = -1.5:0.001:1.5
y_range = -1.5:0.001:1.5
x = repeat(collect(x_range), 1, length(y_range)) 
y = repeat(collect(y_range)', length(x_range), 1)  
z = zeros(size(x))  

# Potential function calculation
r1 = sqrt.((x .+ mu) .^ 2 .+ y .^ 2 .+ z .^ 2)
r2 = sqrt.((x .- (1 - mu)) .^ 2 .+ y .^ 2 .+ z .^ 2)
U = 0.5 .* (x .^ 2 .+ y .^ 2) .+ (1 - mu) ./ r1 .+ mu ./ r2
C = 2 .* U

"""manifolds"""
tf1, tf2, tf3 = 15.0, 25.0, 145.0
sol_sp_L1, sol_sm_L1, sol_up_L1, sol_um_L1 = fun_compute_manifolds(L1, mu, tf1)
sol_sp_L2, sol_sm_L2, sol_up_L2, sol_um_L2 = fun_compute_manifolds(L2, mu, tf2)
sol_sp_L3, sol_sm_L3, sol_up_L3, sol_um_L3 = fun_compute_manifolds(L3, mu, tf3)

# Lagrange Points and Their Corresponding Jacobi Constants
C1 = fun_Jacobi_const(vcat(L1, [0.0, 0.0, 0.0]), mu)
C2 = fun_Jacobi_const(vcat(L2, [0.0, 0.0, 0.0]), mu)
C3 = fun_Jacobi_const(vcat(L3, [0.0, 0.0, 0.0]), mu)

# Generate each figure
for (Li, sol_sp, sol_sm, sol_up, sol_um, Ci, tf, label) in 
    zip([L1, L2, L3], [sol_sp_L1, sol_sp_L2, sol_sp_L3], [sol_sm_L1, sol_sm_L2, sol_sm_L3], 
        [sol_up_L1, sol_up_L2, sol_up_L3], [sol_um_L1, sol_um_L2, sol_um_L3], 
        [C1, C2, C3], [tf1, tf2, tf3], ["L1", "L2", "L3"])

    local fig = Figure(size=(800, 600))
    local ax = Axis(fig[1, 1], aspect = DataAspect(), xlabel = "x [-]", ylabel = "y [-]",
        xlabelsize = 16, ylabelsize = 16,
        xticklabelsize = 16, yticklabelsize = 16,
        xlabelfont = font,
        ylabelfont = font,
        xticklabelfont = font,
        yticklabelfont = font
    )

    p_primary = scatter!(ax, [-mu, 1 - mu], [0, 0], color = :black, markersize = 15)

    p_stable1 = lines!(ax, real(sol_sp[1, :]), real(sol_sp[2, :]), color = :blue, linewidth = 1.5)
    p_stable2 = lines!(ax, real(sol_sm[1, :]), real(sol_sm[2, :]), color = :blue, linewidth = 1.5)

    p_unstable1 = lines!(ax, real(sol_up[1, :]), real(sol_up[2, :]), color = :red, linewidth = 1.5)
    p_unstable2 = lines!(ax, real(sol_um[1, :]), real(sol_um[2, :]), color = :red, linewidth = 1.5)

    p_C = contour!(ax, x_range, y_range, C; levels = [Ci], linewidth = 1.5, color = :black)

    p_lagrange = scatter!(ax, [L1[1], L2[1], L3[1], L4[1], L5[1]],
                              [L1[2], L2[2], L3[2], L4[2], L5[2]],
                              color = :black, marker = :star6, markersize = 15)
    p_Li = scatter!(ax, [Li[1]], [Li[2]], color = :orange, marker = :star6, markersize = 15)

    fig[1, 2] = Legend(fig, [
            p_primary,
            p_stable1,
            p_unstable1,
            p_C,
            p_lagrange,
            p_Li
        ], [
            "Primary bodies",
            "Stable manifolds",
            "Unstable manifolds",
            "Zero-velocity curve",
            "Lagrange points",
            "Target L-point"
        ], 
        labelfont = font)

    # Save the figure
    local save_file_name = "Ex3_manifold_of_$(label)_mu=$(mu)_C=$(Ci)_tf=$(tf)"
    local save_file_name = replace(save_file_name, "." => ",")
    local save_path = joinpath(save_fig_dir, save_file_name * ".png")
    save(save_path, fig)

    # Save trajectory data
    fun_save_data(sol_sp.t, sol_sp.u, save_dat_dir, save_file_name * "_stable_plus")
    fun_save_data(sol_sm.t, sol_sm.u, save_dat_dir, save_file_name * "_stable_minus")

    display(fig)
end

