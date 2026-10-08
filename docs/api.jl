# Shared categorization for the public API pages and the developer appendix.
# The compatibility links in reference.md are a snapshot of the former page;
# keep their anchors when moving an existing entry again.
module DocumentationAPI
using JosephsonCircuits

const GROUPS = (
    circuits = Symbol.([
        "Circuit", "Interface", "Instance", "Ground", "Net", "PortRef", "PinRef",
        "compile", "elaborate", "ElaboratedCircuit", "CompiledCircuit", "symbolicmatrices",
        "numericmatrices", "calccircuitgraph", "IctoLj", "LjtoIc", "@params",
        "ComponentNotSupportedError",
    ]),
    components = Symbol.([
        "Inductor", "Capacitor", "Resistor", "CurrentSource", "VoltageSource", "Port",
        "MatchedTermination", "MutualInductor", "JosephsonJunction", "NonlinearInductor",
        "PolynomialCPR", "FrequencyDependent", "LaplaceResponse", "ScatteringParameters",
        "GaussianChannel", "TransmissionLine", "RationalScattering", "VectorFitting",
        "PassivityEnforcement", "LinearizedScattering", "Passive", "Lossless",
        "ScatteringLimit", "OpenDC", "ShortDC", "ThroughDC", "ScatteringDC",
        "ThermalEquilibrium", "NoiseCovariance", "ConjugateSymmetry", "Native",
        "passivityassessment",
    ]),
    harmonicbalance = Symbol.([
        "hbsolve", "hbnlsolve", "hblinsolve", "hbstability", "HBStabilityResult",
        "Monodromy", "ShiftInvert", "DenseSpectrum", "ContourIntegral", "matchpoles",
        "NonlinearHB", "LinearizedHB", "HB", "designsensitivities", "designjacobian",
        "hbcache", "hbsolve!", "HBCache", "HBReuse", "reset!",
    ]),
    solvers = Symbol.([
        "NewtonKrylov", "Newton", "QuasiNewton", "ExternalSolver", "GMRES", "KrylovJL",
        "Staged", "BlockDiagonal", "FullJacobian", "HarmonicBand", "MeasuredBand",
        "Clusters", "CoupledModes", "CouplingMask", "Automatic", "Floquet", "Always",
        "Probe", "Never", "Backtracking", "KLUfactorization", "LUfactorization",
        "QRfactorization", "CUDSSFactorization", "BlockFactorization",
    ]),
    interop = Symbol.([
        "HBNonlinearProblem", "hbnonlinearproblem", "JacobianOperator", "preconditioner",
        "hbresidual!", "hbjvp!", "hbvjp!", "hbjacobian!", "hbd2F!", "hbd3F!", "hbdFdp!",
        "jacobianprototype", "setdrive!", "drivenresidual!",
    ]),
    transient = Symbol.([
        "TransientSource", "TransientState", "transientproblem", "transientstate",
        "transientsolve", "transientsensitivity", "transientdemodulate",
        "transienttangent", "transientadjoint", "transientinjection", "Trapezoidal",
        "GaussLegendre", "BackwardEuler", "WRspice", "TransientReuse", "TransientSolution",
        "TransientBatchSolution", "TransientStepError",
    ]),
    noise = Symbol.([
        "transientiqplan", "transientiq!", "transientiq", "transientiqvjp!",
        "transientquantumplan", "transientquantum", "transientquantum!",
        "transientquantumvjp!", "transientnoisebaths", "transientnoise", "transientgain",
        "transientquantumdiagnostics", "transientquantumefficiency", "thermaloccupation",
        "effectivetemperature", "noisetemperature", "noisequanta",
    ]),
    networks = Symbol.([
        "connectS", "solveS", "quadraturetransform", "phi0", "Phi0", "speed_of_light",
        "planck_constant", "reduced_planck_constant", "boltzmann_constant",
    ]),
)

in_group(obj, group) = any(name -> isdefined(JosephsonCircuits, name) &&
    getfield(JosephsonCircuits, name) === obj, getproperty(GROUPS, group))
is_internal(obj) = !any(group -> in_group(obj, group), keys(GROUPS))
end
