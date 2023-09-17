using Test, Documenter, JosephsonCircuits

DocMeta.setdocmeta!(JosephsonCircuits, :DocTestSetup, :(using JosephsonCircuits); recursive=true)
makedocs(modules = [JosephsonCircuits], sitename="JosephsonCircuits.jl")

deploydocs(
   repo = "github.com/kpobrien/JosephsonCircuits.jl.git",
)
