using Test, Documenter, JosephsonCircuits

DocMeta.setdocmeta!(JosephsonCircuits, 
	:DocTestSetup, 
	:(using JosephsonCircuits);
	recursive=true)
#makedocs(sitename="My Documentation",modules=[QCE])

doctest(JosephsonCircuits, manual = false)


