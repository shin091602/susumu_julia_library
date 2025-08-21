# =================================
# =================================
#
# Compute orbits
# By: Damennick Henry
# Date: 9/29/22
# 
# Description: Compute qpos around the 9:2 NRHO
#
# Notes: 
#
# =================================
# =================================
using CSV, Tables, LinearAlgebra, SparseArrays, Statistics, CDDLib, Polyhedra, NearestNeighbors, GLPK, TickTock, Dates, DelimitedFiles, Dierckx, JLD2, StaticArrays, Plots, Roots, DifferentialEquations, MAT, Interpolations, BSON

include("cr3bp.jl")
include("po_ms.jl")
include("gmos.jl")
include("genfuncs.jl")
include("qpo_general.jl")
include("pseudoarclength.jl")

function load_parameters()
    # =================================
    # =================================
    #
    # Load Parameters
    # By: Damennick Henry
    # Date: 3/15/22
    # 
    # Description: Load parameters for transfers
    #
    # Notes: 
    #
    # =================================
    # =================================

    p = Dict();

    # SYSTEM PARAMETERS
    # Gravitational parameter
    p["mu"] = 0.012153599037880;
    # Jacobi integral
    p["C"] = 3.046741234093354;
    # State dimension
    p["d"] = 6;
    # State and parameter dimension
    p["dp"] = 6;

    # PO ORBIT PARAMETERS
    # Angles to compute orbit at
    tht0 = LinRange(0, 2*pi, 26); p["tht0"] = 0:step(tht0):tht0[end]-1e-5;
    # Number of invariant curve points
    p["N"] = length(p["tht0"]);
    # Number of multiple shooting points
    p["M"] = 25;
    # Phase condition flag
    p["phsflg"] = true;
    # Parameterization constraints
    p["s"] = [(z, Xd) -> jacobi_integral(z[1:6], p["mu"]) - p["C"]];
    # Parameterization constraint Jacobians
    p["Ds"] = [(z, Xd) -> [DC(z[1:6], p["mu"])' zeros(1, (p["M"]-1)*p["d"] + 1)]];
    # Step size onto manifold
    p["eps"] = 1e-5;
    # Manifold directions
    p["diru"] = [1; 0; 0; 0; 0; 0];
    p["dirs"] = [-1; 0; 0; 0; 0; 0];

    # TORUS COMPUTATION PARAMETERS
    p["NT"] = (25, 25);
    # Dummy dx variable
    p["dx"] = zeros(p["d"]);
    # Invariant curve and initial STM 
    p["uphi0"] = [zeros(p["d"]*p["NT"][2]); repeat(reshape(1.0*Matrix(I, p["dp"], p["dp"]), p["dp"]^2, 1), p["NT"][2])];
    # GMOS eom function
    idx = Vector{Bool}(undef, p["NT"][2]*p["d"] + p["NT"][2]*(p["d"]^2));
    idx .=0;
    p["gEOM"] = (dx, x, q, t) -> eom_gmos!(dx, x, p, t, zeros(42), idx, p);
    # Phase condition index
    p["pconidx"] = 1:2;
    # Linear step on to qpo
    p["K"] = 1e-6;
    # Floquet matrix
    p["B"] = zeros(p["d"]*p["NT"][2], p["d"]*p["NT"][2]);
    # Rotation matrix
    p["Q"] = zeros(p["d"]*p["NT"][2], p["d"]*p["NT"][2])
    # Rotation transformation
    p["R"] = zeros(p["d"]*p["NT"][2], p["d"]*p["NT"][2])
    # Tolerance for eigenvalues for manifold tangent directions
    p["tevtol"] = 1e-3;
    # Linear step size on to manifold
    p["eps"] = 1e-5;
    # Center tanget switch
    p["ctswitch"] = true;
    # Number of time points along manifold
    p["tp"] = 1000;
    
    
    # VECTOR FIELD FUNCTIONS
    # Define the EOM functions
    p["f"] = (dx, x, q, t) -> f_CR3BP!(dx, x, p["mu"]);
    p["df"] = (df, x, q, t) -> Df_CR3BP!(df, x, p["mu"], p["d"]);
    p["eom"] = (dx, x, q, t) -> eom!(dx, x, p, t, zeros(6, 6), p["f"], p["df"], p["d"]);
    p["f_ivc!"] = (dx, x, q, t) -> vec(vcat([p["f"](@view(dx[i*p["d"] - (p["d"]-1):i*p["d"]]), x[i*p["d"] - (p["d"]-1):i*p["d"]], q, t) for i in 1:p["NT"][2]]...));

    return p
end

function main()
    # =================================
    # =================================
    #
    # Main
    # By: Damennick Henry
    # Date: 3/15/22
    # 
    # Description: Main function for computing transfers
    #
    # Notes: 
    #
    # =================================
    # =================================
    gr()
    # INITIALIZATION
    # Load parameters
    p = load_parameters();

    # COMPUTE PERIODIC ORBIT
    zpo = matread(joinpath(@__DIR__, "NRHO92.mat"));
    zpo = [zpo["x0"]; zpo["T"]];
    # Initialize for multiple shooting
    (z0, UdA) = po_ms_initialization(zpo, p);
    # Refine orbit
    zpoA = po_refine(z0, UdA, p);
 
    zpoB = copy(zpoA); UdB = deepcopy(UdA)

    # COMPUTE QPO
    # Angles to compute orbit at
    tht0 = LinRange(0, 2*pi, 26); p["tht0"] = 0:step(tht0):tht0[end]-1e-5;
    # Number of invariant curve points
    p["N"] = 25;
    # Fourier matrices
    (p["Fr"], p["IFr"], p["DFr"]) = fourier_matrix(p["d"], p["NT"][2]);
    (p["Fr0"], p["IFr0"], p["DFr0"]) = fourier_matrix(p["d"], p["NT"][1]);
    (p["Fr1"], p["IFr1"], p["DFr1"]) = fourier_matrix(p["d"], p["NT"][2]);
    # Parameterization constraints (fix the stroboscopic time)
    p["s"] = [(z, Ud) -> z[end-3] - zpo[end]];
    # Parameterization constraint Jacobians (fix the stroboscopic time)
    p["Ds"] = [(z, Ud) -> [zeros(1, p["N"]*p["M"]*p["d"]) 1 0 0 0]];
    # Error vector and Jacobian functions
    F(z, Ud) = F_gmos(z, Ud, p);
    DF(z, Ud) = DF_gmos(z, Ud, p);
    Udf!(z, Ud) = Ud_finalization_gmos!(z, Ud, p);
    # Continuation parameters
    pac_p = pac_parameters(100, 1e-12, 20, 5, 1e-5, 1:p["d"]*p["N"]*p["M"]+4, true, p["N"]*p["M"], p["N"]*p["M"]*p["d"], false, false, false);
    # Compute the center tangent direction and torus
    (z0, phi0, Ud0) = center_tangent_qpo([zpoB[1:6]; zpoB[end]], p);
    Z = pa_cont(z0, 1e-6, Ud0, phi0, F, DF, Udf!, pac_p);
    zB  = Z;

    # SAVE DATA
    # Save orbit in MAT file
    matwrite(joinpath(@__DIR__, "orbits.mat"), Dict(
        "UB" => zB[end].z
    ))
    

    UdA = nothing;
    UdB = nothing;
    GC.gc();
end

main();
