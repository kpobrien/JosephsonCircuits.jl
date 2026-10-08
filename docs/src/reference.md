# API reference

Choose a topic for its public interfaces. For units and array conventions,
see [Conventions](conventions.md); for a runnable calculation, start with
the [quickstart](quickstart.md). Public names which are not exported, such
as `JosephsonCircuits.@params`, `JosephsonCircuits.reset!` and the network
library, are written qualified with the module's name and listed with the
related exported functions.

- [Circuit construction and physical constants](api/circuits.md)
- [Components and scattering models](api/components.md)
- [Harmonic-balance analyses](api/harmonicbalance.md)
- [Nonlinear and linear solver options](api/solvers.md)
- [External solver interface](api/interop.md)
- [Transient analyses and responses](api/transient.md)
- [Temporal measurements and noise](api/noise.md)
- [Network parameters](api/networkparameters.md)
- [Closed-form networks](api/networkmodels.md)
- [Network connections](api/connections.md)

The [developer appendix](api/internals.md) contains the remaining internal
and lower-level interfaces, separate from the public reference.
