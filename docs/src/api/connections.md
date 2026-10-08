# Network connections

Join scattering networks at their ports, with their noise: `connectS`
connects them one pair of ports at a time, `solveS` solves the whole
connection at once, and `cascadeS`, `intraconnectS` and `interconnectS`
are the operations they are built from. The noise covariances are
symmetrized and counted in quanta, the vacuum being half a photon: a
passive network given without one emits `(I - S*S')/2`. See
[conventions](../conventions.md#Noise-normalization-and-temperature).

See [network parameters](networkparameters.md) and
[closed-form networks](networkmodels.md), or the
[API overview](../reference.md) for other topics.

```@autodocs
Modules = [JosephsonCircuits]
Filter = obj -> Main.DocumentationAPI.in_group(obj, :connections)
```
