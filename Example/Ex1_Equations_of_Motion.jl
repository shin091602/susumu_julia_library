"""
Ex1_Equations_of_Motion.jl
"""

using LinearAlgebra
using DifferentialEquations
using GLMakie 

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

# Start of script
# closeall() # Close all figures
# Start timer
t0 = time()  
# font setting
font = "Times New Roman"
# Save Figure directory
save_fig_dir = "Figure/"
# Save dat directory
save_dat_dir = "Output_data/"

"""Initial conditions and constants"""
# Mass of Sun [kg]
M_S = 1.9885e30
# Mass of Earth [kg]  
M1 = 5.9724e24
# Mass of Moon [kg]
M2 = 7.346e22
# Gravitational constant [m^3/(kg*s^2)]
G = 6.673e-11

"""Earth-Moon CR3BP"""
# Characteristic length of Earth-Moon CR3BP [m] 
chara_length_CR3BP = 3.8440e8 
# Characteristic mass of Earth-Moon CR3BP [kg] 
chara_mass_CR3BP = M1 + M2
# Characteristic time of Earth-Moon CR3BP [s] 
chara_time_CR3BP = sqrt(chara_length_CR3BP^3 / (G * chara_mass_CR3BP))
# Mass ratio of Earth-Moon CR3BP [no unit]
mu_EM = M2 / (M1 + M2)
# Mean-motion of Earth-Moon CR3BP [1/s]
N_CR3BP = sqrt(G * chara_mass_CR3BP / chara_length_CR3BP^3)

"""Earth-Moon ER3BP"""
# Semimajor axis of Earth-Moon ER3BP [m]
a = 3.8440e8
# Mean-motion of Earth-Moon ER3BP [1/s]
N_ER3BP = sqrt(G * (M1 + M2) / (a^3))
# Eccentricity of Earth-Moon ER3BP [no unit]
e = 0.0549 

"""Earth-Moon-Sun BCR4BP"""
# Characteristic length of the Sun-B1 [m] 
chara_length_SB1 = 1.4960e11
# Characteristic mass of the Sun-B1 [kg] 
chara_mass_SB1 = M1 + M2 + M_S 
# Characteristic time of the Sun-B1 [s] 
chara_time_SB1 = sqrt(chara_length_SB1^3 / (G * chara_mass_SB1))
# Nondimensional mass of Sun [s]
m_S = M_S / (M1 + M2)
# Nondimensional Sun orbit radius [m]
a_S = chara_length_SB1 / chara_length_CR3BP
# Nondimensional Sun angular velocity [1/s]
omega_S = sqrt((1 + m_S) / (a_S^3)) - 1
# Initial phase angle of the Sun-B1 [rad]
theta_S0 = π 
# Mass ratio of the Sun-B1 [no unit]
mu_SB1 = (M1 + M2) / (M1 + M2 + M_S)
# Nondimensional distance of Earth-Moon-Sun BCR4BP [no unit]
a_EM = chara_length_CR3BP / chara_length_SB1
# nondimensional Moon angular velocity [no unit]
omega_M = chara_time_SB1 / chara_time_CR3BP - 1
# Initial phase angle of Earth-Moon [rad]
theta_M0 = 0

"""L2 Lyapunov orbit in Earth-Moon CR3BP"""
# Initial Condition 
t_n_CR3BP = 3.467622949281189
x_n_CR3BP = [1.102866098413080; 0.0; 0.0; 0.0; 0.259155907029058; 0.0]

t_CR3BP = t_n_CR3BP
x0_CR3BP = x_n_CR3BP

# Specification Parameter
parameter1 = (mu_EM)
# Specification Integration Time 
tspan1 = (0.0, t_CR3BP)
# Setup ODE Problem
prob1 = ODEProblem(fun_cr3bp!, x0_CR3BP, tspan1, parameter1)
# ODE Solver Execution
sol1 = solve(prob1, Vern7(), abstol=1e-14, reltol=1e-14)

save_file_name1 = "Ex1_CR3BP_x0=$(x0_CR3BP[1])_vy0=$(x0_CR3BP[5])_t=$(t_CR3BP)"
save_file_name1 = replace(save_file_name1, "." => ",")

fun_make_fig_2D_orbit(sol1[1, :], sol1[2, :], "x", "y", 2, 20, save_fig_dir, save_file_name1, mu_EM)
fun_save_data(sol1.t, sol1.u, save_dat_dir, save_file_name1)

"""trajectory in Earth-Moon ER3BP"""
# Orbital Period of ER3BP
t_ER3BP = t_n_CR3BP * chara_time_CR3BP * N_ER3BP
# transfomation from CR3BP reference frame to the inertial frame
X0_CR3BP = [x_n_CR3BP[1:3] .* chara_length_CR3BP; x_n_CR3BP[4:6] .* chara_length_CR3BP ./ chara_time_CR3BP]
C_CR3BP = [cos(0.0) -sin(0.0) 0.0;
           sin(0.0)  cos(0.0) 0.0;
                0.0       0.0 1.0]
dtheta_dt_CR3BP = N_CR3BP
# constructs a transformation matrix of CR3BP
rotating_matrix_CR3BP = fun_rotating_to_inertial_matrix(C_CR3BP, dtheta_dt_CR3BP)
# Convert the initial state vector X0_CR3BP expressed in the rotating frame of CR3BP to the inertial frame
X0_inertial = rotating_matrix_CR3BP * X0_CR3BP

# Transformation from inertial frame to ER3BP
C_ER3BP = [cos(0.0) -sin(0.0) 0.0;
           sin(0.0)  cos(0.0) 0.0;
                0.0       0.0 1.0]
# Calculation of angular velocity of ER3BP
dtheta_dt_ER3BP = sqrt(G * (M1 + M2) * (1 + e * cos(0.0))^4 / (a * (1 - e^2))^3)
# constructs a transformation matrix of ER3BP
rotating_matrix_ER3BP = fun_rotating_to_inertial_matrix(C_ER3BP, dtheta_dt_ER3BP)
# Convert the initial state vector X0_inertial expressed in the inertial frame to the rotating frame of ER3BP
X0_ER3BP = rotating_matrix_ER3BP \ X0_inertial # (inv(rotating_matrix_ER3BP)*X0_inertialと同じ意味)
x0_ER3BP = [X0_ER3BP[1:3] ./ (a * (1 - e)); X0_ER3BP[4:6] ./ (a * (1 - e)) ./ N_ER3BP] # x0とX0で座標系が異なるので注意。

# Specification Parameter 
parameter2 = (mu_EM, e)
# Specification Integration Time 
tspan2 = (0.0, t_ER3BP)
# Setup ODE Problem 
prob2 = ODEProblem(fun_er3bp!, x0_ER3BP, tspan2, parameter2)
# ODE Solver Execution
sol2 = solve(prob2, Vern7(), abstol=1e-14, reltol=1e-14)

# Save figure
save_file_name2 = "Ex1_ER3BP_x0=$(x0_ER3BP[1])_vy0=$(x0_ER3BP[5])_t=$(t_ER3BP)"
save_file_name2 = replace(save_file_name2, "." => ",")

fun_make_fig_2D_orbit(sol2[1, :], sol2[2, :], "x", "y", 2, 20, save_fig_dir, save_file_name2, mu_EM)
fun_save_data(sol2.t, sol2.u, save_dat_dir, save_file_name2)

"""trajectory in Earth-Moon_Sun BCR4BP in the Earth-Moon rotating frame"""
s# Define the integration time for the BCR4BP
t_BCR4BP = t_n_CR3BP * chara_time_CR3BP * N_CR3BP

# Specification Parameter for BCR4BP in the Earth-Moon rotating frame
parameter3 = (mu_EM, m_S, a_S, omega_S, theta_S0)

# Specification integration time 
tspan3 = (0.0, t_BCR4BP)
# Setup ODE Problem 
prob3 = ODEProblem(fun_bcr4bp_EMS!, x0_CR3BP, tspan3, parameter3)
# ODE Solver Execution
sol3 = solve(prob3, Vern7(), abstol=1e-14, reltol=1e-14)

# Save the figure
save_file_name3 = "Ex1_BCR4BP_EMS_x0=$(x0_ER3BP[1])_vy0=$(x0_ER3BP[5])_t=$(t_BCR4BP)"
save_file_name3 = replace(save_file_name3, "." => ",")

fun_make_fig_2D_orbit(sol3[1, :], sol3[2, :], "x", "y", 2, 20, save_fig_dir, save_file_name3, mu_EM)
fun_save_data(sol3.t, sol3.u, save_dat_dir, save_file_name3)

"""trajectory in Earth-Moon_Sun BCR4BP in the Sun-B1 rotating frame"""
# Define angle range
theta = range(0, 2 * pi, length=1001)

# Compute Moon and Earth trajectories
x_moon = 1 - mu_SB1 .+ a_EM * (1 - mu_EM) .* cos.(theta)
y_moon = a_EM * (1 - mu_EM) .* sin.(theta)
x_Earth = 1 - mu_SB1 .+ a_EM * mu_EM .* cos.(theta)
y_Earth = a_EM * mu_EM .* sin.(theta)

# Time scaling transformation
t_bcr4bp_SB1 = t_CR3BP * chara_time_CR3BP / chara_time_SB1

# Transformation from the inertial frame to the Sun-B1 rotating frame
C_SB1 = [cos(0) -sin(0) 0;
         sin(0)  cos(0) 0;
              0       0 1]

dtheta_dt_SB1 = 1 / chara_time_SB1
rotating_matrix_SB1 = fun_rotating_to_inertial_matrix(C_SB1, dtheta_dt_SB1)

# Convert initial conditions from the inertial frame to the Sun-B1 rotating frame
X0_bcr4bp_SB1 = rotating_matrix_SB1 \ X0_inertial
x0_bcr4bp_SB1 = vcat(X0_bcr4bp_SB1[1:3] ./ chara_length_SB1, X0_bcr4bp_SB1[4:6] ./ chara_length_SB1 .* chara_time_SB1)

# Apply parallel translation
x0_bcr4bp_SB1[1] += (1 - mu_SB1)

# Define time span
tspan4 = (0.0, t_bcr4bp_SB1)

# Specification Parameter for BCR4BP in the Sun-B1 rotating frame
parameter4 = (mu_SB1, mu_EM, a_S, a_EM, omega_M, theta_M0)
parameter4 = (mu_SB1, mu_EM, 1, a_EM, omega_M, theta_M0)

# Solve the ODE
prob4 = ODEProblem(fun_bcr4bp_SB1!, x0_bcr4bp_SB1, tspan4, parameter4)
sol4 = solve(prob4, Vern7(); reltol=1e-9, abstol=1e-9)

# Save the figure
save_file_name4 = "Ex1_BCR4BP_SB1_x0=$(x0_ER3BP[1])_vy0=$(x0_ER3BP[5])_t=$(t_BCR4BP)"
save_file_name4 = replace(save_file_name4, "." => ",")

fig4 = Figure(size = (800, 600))

ax4 = Axis(fig4[1, 1],
    xlabel = "x [-]", ylabel = "y [-]",
    aspect = DataAspect(),
    xlabelsize = 20, ylabelsize = 20,
    xticklabelsize = 20, yticklabelsize = 20,
    xlabelfont = font,
    ylabelfont = font,
    xticklabelfont = font,
    yticklabelfont = font
)

# Earth and Moon trajectories (if defined)
lines!(ax4, x_moon, y_moon, color = "#EDB120", label = "Moon Orbit")
lines!(ax4, x_Earth, y_Earth, color = "#7E2F8E", label = "Earth Orbit")
# Initial spacecraft position
scatter!(ax4, [sol4[1, 1]], [sol4[2, 1]], color = "#0072BD", markersize = 10, label = "Initial Position")
# Spacecraft orbit
lines!(ax4, sol4[1, :], sol4[2, :], color = "#0072BD", linewidth = 2, label = "Orbit")
# Legend
fig4[1, 2] = Legend(fig4, ax4, "Legend", orientation = :vertical, labelfont = font)
# Save figure
save(joinpath(save_fig_dir, save_file_name4 * ".png"), fig4)
fig4

fun_save_data(sol4.t, sol4.u, save_dat_dir, save_file_name4)