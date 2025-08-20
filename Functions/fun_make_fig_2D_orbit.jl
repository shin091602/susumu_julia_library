function fun_make_fig_2D_orbit(data_x, data_y, xlabel_name, ylabel_name, orbit_line_width, label_font_size, save_dir, save_file_name, mu_EM)
    # font setting
    font = "Times New Roman"
    
    fig = Figure(size = (800, 600))
    
    # Create an Axis for the Plot
    ax = Axis(fig[1, 1], 
    xlabel = xlabel_name, ylabel = ylabel_name, 
    aspect = DataAspect(), xlabelsize = label_font_size, ylabelsize = label_font_size,
    xticklabelsize = label_font_size, yticklabelsize = label_font_size,
    xlabelfont = font,
    ylabelfont = font,
    xticklabelfont = font,
    yticklabelfont = font
    )
    
    # Plot the Moon's position
    scatter!(ax, [1 - mu_EM], [0.0], color = "#EDB120", markersize = 20, label = "Moon")
    
    # Plot the initial position of the spacecraft
    scatter!(ax, [data_x[1]], [data_y[1]], color = "#0072BD", markersize = 10, label = "Initial Position")
    
    # Plot the trajectory of the spacecraft
    lines!(ax, data_x, data_y, color = "#0072BD", linewidth = orbit_line_width, label = "Orbit")
    
    # Add a Legend
    legend = Legend(fig, ax, "Legend", orientation = :vertical, labelfont = font)
    fig[1, 2] = legend
    
    # Save figure
    save_path = joinpath(save_dir, save_file_name * ".png")
    save(save_path, fig)
    
    # Display the Figure
    fig
end
