using Pkg
Pkg.activate(@__DIR__)
Pkg.add("JuliaFormatter")
using JuliaFormatter

format("src")
println("Done.")
