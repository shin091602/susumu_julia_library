"""
Ex2_Zero_velocity_curve.jl
"""

using GLMakie
using Polynomials
using CairoMakie 

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
# Save Figure directory
save_fig_dir = "Figure/"

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

# Plot
fig = Figure(size=(600, 600))

# Create an Axis for the Plot
ax = Axis(fig[1, 1], aspect = DataAspect(), xlabel = "x [-]", ylabel = "y [-]",
    xlabelsize = 16, ylabelsize = 16,
    xticklabelsize = 16, yticklabelsize = 16,
    xlabelfont = font,
    ylabelfont = font,
    xticklabelfont = font,
    yticklabelfont = font
)

# Define contour levels
Clevels = [3.000, 3.020, 3.040, 3.060, 3.080, 3.100, 3.120, 3.140, 3.160, 3.180, 3.200, 3.220, 3.240, 3.260, 3.280, 3.300]

# Get a colormap with the same number of levels
colors = cgrad(:jet, length(Clevels), categorical=true)  # Get distinct colors for each level

# Manually plot each contour level with a different color
for (i, level) in enumerate(Clevels)
    GLMakie.contour!(ax, x_range, y_range, C; levels=[level], linewidth=1.5, color=colors[i])
end

# Plot primary celestial bodies
GLMakie.scatter!(ax, [-mu, 1 - mu], [0, 0], color=:black, markersize=10)

# Plot Lagrange points
GLMakie.scatter!(ax, [L1[1], L2[1], L3[1], L4[1], L5[1]], [L1[2], L2[2], L3[2], L4[2], L5[2]], 
         color=:black, marker=:star5, markersize=10)

# Colorbar (ensures colors are mapped correctly)
Colorbar(fig[1, 2], colormap=:jet, limits=(minimum(Clevels), maximum(Clevels)), ticks=Clevels, 
    label = "Jacobi constant [-]", 
    labelfont = font, 
    ticklabelfont = font
)

display(fig)

# Save figure
save_file_name = "Ex2_Zero_Velocity_Curve_mu$(mu)"
save_path = joinpath(save_fig_dir, save_file_name * ".png")
save(save_path, fig)

