# Circuit construction

Define topology and bind component values. Start with the [circuit guide](../circuits.md).

See the [API overview](../reference.md) for other analyses and
[conventions](../conventions.md) for units and array ordering.

```@autodocs
Modules = [JosephsonCircuits, JosephsonCircuits.CircuitValues]
Filter = obj -> Main.DocumentationAPI.in_group(obj, :circuits)
```

## Physical constants

The flux quanta follow from the exact SI values of the Planck constant
and the elementary charge.

```@autodocs
Modules = [JosephsonCircuits]
Filter = obj -> Main.DocumentationAPI.in_group(obj, :constants)
```
