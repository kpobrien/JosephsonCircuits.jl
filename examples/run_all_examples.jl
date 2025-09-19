# examples/run_all_examples.jl
using Pkg
Pkg.activate(@__DIR__)  # Use the examples environment
Pkg.instantiate()  # Install dependencies if needed

using JosephsonCircuits
using CairoMakie

println("Running all NL element examples...")
println("="^50)

examples = [
    "01_simple_jpa_example.jl",
    "02_jtwpa_example.jl", 
    "03_snailpa_example.jl",
    "04_flux_jtwpa_example.jl",
    "05_floquet_jtwpa_loss_example.jl"
]

for example in examples
    println("\nRunning $example...")
    include(example)
    println("\n" * "="^50)
end

println("\nAll examples completed!")
println("Check the generated .png files for comparison plots.")