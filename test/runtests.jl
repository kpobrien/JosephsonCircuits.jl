using Test, Documenter, JosephsonCircuits



DocMeta.setdocmeta!(JosephsonCircuits, :DocTestSetup, :(using JosephsonCircuits,SparseArrays,BandedMatrices,LinearAlgebra); recursive=true)
#makedocs(sitename="My Documentation",modules=[QCE])

# one issue with the FFT tests is that they error if the FFTW provider
# is the default fftw instead of mkl. This can be changed to MKL with the
# command FFTW.set_provider!("mkl")

doctest(JosephsonCircuits; manual = false)

