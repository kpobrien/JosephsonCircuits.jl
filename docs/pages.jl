# Shared by the production build and the Node-free documentation check.
docpages = [
    "Start here" => [
        "Home" => "index.md",
        "First calculation" => "quickstart.md",
        "Units and conventions" => "conventions.md",
    ],
    "Build a circuit" => [
        "Components and connections" => "circuits.md",
        "Scattering blocks" => "scattering.md",
        "Migration guide" => "migration.md",
    ],
    "Run an analysis" => [
        "Harmonic balance" => "harmonicbalance.md",
        "Noise at the ports" => "portnoise.md",
        "Transient simulation" => "transient.md",
        "Quantum noise in time" => "transientnoise.md",
    ],
    "Examples" => [
        "Choose an example" => "examples.md",
        "Current-pumped JPA" => "recipes/jpa.md",
        "Environmental impedance" => "recipes/environment.md",
        "Double-pumped JPA" => "recipes/double-pump.md",
        "Flux-pumped JPA" => "recipes/flux-pump.md",
        "SNAIL amplifier" => "recipes/snail.md",
        "Traveling-wave amplifier" => "recipes/traveling-wave.md",
        "Floquet JTWPA" => "recipes/floquet.md",
        "Impedance-engineered JPA" => "recipes/lesa.md",
        "Design sensitivities" => "recipes/sensitivities.md",
        "Scattering sensitivities" => "recipes/scattering-sensitivities.md",
        "Direct current" => "recipes/dc.md",
        "Pulsed Josephson line" => "recipes/transient-line.md",
        "Transient WRspice comparison" => "recipes/transient-wrspice.md",
    ],
    "Advanced use" => [
        "Performance and reuse" => "performance.md",
        "Using other solvers" => "interop.md",
    ],
    "Theory and development" => [
        "Harmonic-balance theory" => "harmonicbalancetheory.md",
        "Transient theory" => "transienttheory.md",
        "Implementation notes" => "implementation.md",
        "Numerical references" => "numerical-references.md",
    ],
    "API reference" => "reference.md",
]
