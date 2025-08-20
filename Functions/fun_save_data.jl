using DataFrames, CSV

function fun_save_data(data_t, data_u, save_dir, save_file_name)
    # Convert sol.u (Vector of Vectors) to Matrix if needed
    if eltype(data_u) <: AbstractVector  # e.g., Vector{Vector{Float64}}
        data_u_mat = hcat(data_u...)'  # Convert to Matrix: (N × dim)
    elseif isa(data_u, AbstractMatrix)
        data_u_mat = data_u
    else
        error("Unsupported data_u type: must be Vector of Vectors or Matrix")
    end

    # Create column names
    if size(data_u_mat, 2) == 6
        col_names = ["t", "x", "y", "z", "dx", "dy", "dz"]
    else
        col_names = ["t"; ["x[$i]" for i in 1:size(data_u_mat, 2)]]
    end

    # Combine time and state into a DataFrame
    df = DataFrame(hcat(data_t, data_u_mat), Symbol.(col_names))

    # Ensure directory exists
    dir_path = dirname(save_file_name)
    isdir(dir_path) || mkpath(dir_path)

    # Sanitize file name (optional)
    safe_name = replace(save_file_name, "," => "_")

    # Write to .dat file (CSV format)
    CSV.write(joinpath(save_dir*save_file_name*".dat"), df)
end

