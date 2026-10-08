# Shared categorization for the public API pages and the developer appendix.
module DocumentationAPI
using JosephsonCircuits
const JC = JosephsonCircuits

const GROUPS = (
    circuits = Symbol.([
        "Circuit", "Interface", "Ground", "Net", "PortRef", "PinRef",
        "compile", "CompiledCircuit", "symbolicmatrices", "numericmatrices",
        "IctoLj", "LjtoIc", "@params",
        "ComponentNotSupportedError",
    ]),
    constants = Symbol.([
        "phi0", "Phi0", "elementary_charge", "speed_of_light", "planck_constant",
        "reduced_planck_constant", "boltzmann_constant",
    ]),
    components = Symbol.([
        "Inductor", "Capacitor", "Resistor", "CurrentSource", "Port",
        "MatchedTermination", "MutualInductor", "JosephsonJunction", "NonlinearInductor",
        "PolynomialCPR", "FrequencyDependent", "LaplaceResponse", "ScatteringParameters",
        "TransmissionLine", "RationalScattering", "VectorFitting",
        "PassivityEnforcement", "LinearizedScattering", "Passive", "Lossless",
        "ScatteringLimit", "OpenDC", "ShortDC", "ThroughDC", "ScatteringDC",
        "ThermalEquilibrium", "NoiseCovariance", "ConjugateSymmetry", "Native",
        "passivityassessment",
    ]),
    harmonicbalance = Symbol.([
        "hbsolve", "hbnlsolve", "hblinsolve", "hbstability", "HBStabilityResult",
        "Monodromy", "ShiftInvert", "DenseSpectrum", "ContourIntegral", "matchpoles",
        "NonlinearHB", "LinearizedHB", "HB", "designsensitivities", "designjacobian",
        "hbcache", "hbsolve!", "HBCache", "reset!",
    ]),
    solvers = Symbol.([
        "NewtonKrylov", "Newton", "QuasiNewton", "ExternalSolver", "GMRES", "KrylovJL",
        "Staged", "MeasuredBand", "Automatic", "Always", "Probe", "Never",
        "Backtracking", "KLUfactorization", "LUfactorization", "QRfactorization",
        "CUDSSFactorization", "BlockFactorization",
    ]),
    interop = Symbol.([
        "HBNonlinearProblem", "hbnonlinearproblem", "JacobianOperator", "preconditioner",
        "hbresidual!", "hbjvp!", "hbvjp!", "hbjacobian!", "hbd2F!", "hbd3F!", "hbdFdp!",
        "setdrive!", "drivenresidual!",
    ]),
    transient = Symbol.([
        "TransientSource", "TransientState", "transientproblem", "transientstate",
        "transientsolve", "transientsensitivity", "transientdemodulate",
        "transienttangent", "transientadjoint", "Trapezoidal", "GaussLegendre",
        "BackwardEuler", "WRspice", "TransientReuse", "TransientSolution",
        "TransientBatchSolution", "TransientStepError",
    ]),
    noise = Symbol.([
        "transientiqplan", "transientiq!", "transientiq", "transientiqvjp!",
        "transientquantumplan", "transientquantum", "transientquantum!",
        "transientquantumvjp!", "transientnoisebaths", "transientnoise", "transientgain",
        "transientquantumefficiency", "thermaloccupation", "effectivetemperature",
        "noisetemperature", "noisequanta",
    ]),
)

# The network library is grouped by the source files that document it: a
# group holds every public name a docstring of one of its files documents,
# so that the module's `public` declaration decides what they show and
# their helpers stay internal.
const FILEGROUPS = (
    networkparameters = ["networks/parameters.jl"],
    networkmodels = ["networks/networks.jl"],
    connections = ["networks/connections.jl"],
)

# a name of the package a user may rely on: exported, or declared public,
# which Julia records from 1.11
ispublicname(name::Symbol) = Base.isexported(JC, name) ||
    (isdefined(Base, :ispublic) && Base.ispublic(JC, name))

# the bindings of the package which carry a docstring, with the source
# files of their docstrings
const DOCUMENTED = [(binding.var, Set(String(get(d.data, :path, ""))
    for d in values(multidoc.docs)))
    for (binding, multidoc) in Base.Docs.meta(JC) if binding.mod === JC]

const FILEGROUPOBJECTS = map(FILEGROUPS) do files
    [getfield(JC, name) for (name, paths) in DOCUMENTED
        if isdefined(JC, name) && ispublicname(name) &&
            any(p -> any(f -> endswith(p, f), files), paths)]
end

function in_group(obj, group)
    if haskey(FILEGROUPOBJECTS, group)
        return any(o -> o === obj, getproperty(FILEGROUPOBJECTS, group))
    end
    return any(name -> isdefined(JC, name) && getfield(JC, name) === obj,
        getproperty(GROUPS, group))
end
is_internal(obj) = !any(group -> in_group(obj, group),
    (keys(GROUPS)..., keys(FILEGROUPS)...))

# Every name of a public page is a public name of the package, and every
# documented public name is on a public page.
if isdefined(Base, :ispublic)
    let grouped = reduce(vcat, collect(GROUPS))
        notpublic = filter(name -> !ispublicname(name), grouped)
        isempty(notpublic) || error("names on the public API pages that are neither exported nor declared public: $(notpublic)")
        ungrouped = [name for (name, _) in DOCUMENTED
            if name !== :JosephsonCircuits && isdefined(JC, name) &&
                ispublicname(name) && is_internal(getfield(JC, name))]
        isempty(ungrouped) || error("documented public names on no public API page: $(ungrouped)")
    end
end
end
