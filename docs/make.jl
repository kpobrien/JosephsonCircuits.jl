using Test, Documenter, JosephsonCircuits

DocMeta.setdocmeta!(JosephsonCircuits, :DocTestSetup, :(using JosephsonCircuits); recursive=true)
#makedocs(sitename="My Documentation",modules=[JosephsonCircuits])
#doctest(QCE; manual = false)
makedocs(sitename="JosephsonCircuits")

#using Documenter
#using JosephsonCircuits

#DocMeta.setdocmeta!(MyPackage, :DocTestSetup, :(using JosephsonCircuits); recursive=true)
#makedocs(modules=[MyPackage], ...)

#makedocs(sitename="My Documentation")
