# Closed-form networks

The parameters of standard networks in closed form: transmission lines,
series and shunt elements, T and Π sections, attenuators, couplers,
circulators, terminations and coupled lines, with the conversions between
the even and odd, Maxwell and mutual descriptions of coupled lines. A
network whose element values, an impedance, an admittance or an
electrical length, are arrays over frequency gives one matrix per
frequency, along the dimensions after the first two.

See [network parameters](networkparameters.md) to convert them and
[connections](connections.md) to join them, or the
[API overview](../reference.md) for other topics.

```@autodocs
Modules = [JosephsonCircuits]
Filter = obj -> Main.DocumentationAPI.in_group(obj, :networkmodels)
```
