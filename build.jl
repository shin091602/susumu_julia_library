# build.jl
using Pkg

packages = [
    "LinearAlgebra",
    "DifferentialEquations",
    "GLMakie",
    "CairoMakie",
    "Polynomials",
    "Printf",
    "Colors",
    "Interpolations",
    "ColorSchemes"
]

for pkg in packages
    println("Adding package: $pkg")
    Pkg.add(pkg)
end

println("Package installation complete.")
