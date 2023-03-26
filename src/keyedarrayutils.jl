

"""
    Stokeyed(S, outputmodes, outputportnumbers, inputmodes, 
        inputportnumbers, w)

Convert a scattering parameter array to a keyed array. 

# Examples
```jldoctest
julia> JosephsonCircuits.Stokeyed([11 12;21 22;;;],[(0,)],[1,2],[(0,)],[1,2],[1.0])
5-dimensional KeyedArray(NamedDimsArray(...)) with keys:
↓   outputmode ∈ 1-element Vector{Tuple{Int64}}
→   outputport ∈ 2-element Vector{Int64}
◪   inputmode ∈ 1-element Vector{Tuple{Int64}}
▨   inputport ∈ 2-element Vector{Int64}
▨   freqindex ∈ 1-element UnitRange{Int64}
And data, 1×2×1×2×1 Array{Int64, 5}:
[:, :, 1, 1, 1] ~ (:, :, (0,), 1, 1):
          (1)  (2)
   (0,)    11   21

[:, :, 1, 2, 1] ~ (:, :, (0,), 2, 1):
          (1)  (2)
   (0,)    12   22
```
"""
function Stokeyed(S, outputmodes, outputportnumbers, inputmodes,
    inputportnumbers, w)
    Nfrequencies = length(w)
    
    return AxisKeys.KeyedArray(
        reshape(S, length(outputmodes), length(outputportnumbers),
            length(inputmodes), length(inputportnumbers), Nfrequencies),
        outputmode = outputmodes,
        outputport = outputportnumbers,
        inputmode = inputmodes,
        inputport = inputportnumbers,
        freqindex=1:Nfrequencies,
    )
end

"""
    Stokeyed(S, outputmodes, outputportnumbers, inputmodes, inputportnumbers)

Convert a scattering parameter array to a keyed array. 

# Examples
```jldoctest
julia> JosephsonCircuits.Stokeyed([11 12;21 22],[(0,)],[1,2],[(0,)],[1,2])
4-dimensional KeyedArray(NamedDimsArray(...)) with keys:
↓   outputmode ∈ 1-element Vector{Tuple{Int64}}
→   outputport ∈ 2-element Vector{Int64}
◪   inputmode ∈ 1-element Vector{Tuple{Int64}}
▨   inputport ∈ 2-element Vector{Int64}
And data, 1×2×1×2 Array{Int64, 4}:
[:, :, 1, 1] ~ (:, :, (0,), 1):
          (1)  (2)
   (0,)    11   21

[:, :, 1, 2] ~ (:, :, (0,), 2):
          (1)  (2)
   (0,)    12   22
```
"""
function Stokeyed(S, outputmodes, outputportnumbers, inputmodes,
    inputportnumbers)
    
    return AxisKeys.KeyedArray(
        reshape(S, length(outputmodes), length(outputportnumbers),
            length(inputmodes), length(inputportnumbers)),
        outputmode = outputmodes,
        outputport = outputportnumbers,
        inputmode = inputmodes,
        inputport = inputportnumbers,
    )
end

"""
    CMtokeyed(CM, outputmodes, outputportnumbers, w)

Convert a commutation relations array to a keyed array.

# Examples
```jldoctest
julia> JosephsonCircuits.CMtokeyed([1 2;3 4;;;],[(0,)],[1,2],[1.0,1.1])
3-dimensional KeyedArray(NamedDimsArray(...)) with keys:
↓   outputmode ∈ 1-element Vector{Tuple{Int64}}
→   outputport ∈ 2-element Vector{Int64}
◪   freqindex ∈ 2-element UnitRange{Int64}
And data, 1×2×2 Array{Int64, 3}:
[:, :, 1] ~ (:, :, 1):
          (1)  (2)
   (0,)     1    3

[:, :, 2] ~ (:, :, 2):
          (1)  (2)
   (0,)     2    4
```
"""
function CMtokeyed(CM, outputmodes, outputportnumbers, w)
    Nfrequencies = length(w)
    
    return AxisKeys.KeyedArray(
        reshape(CM, length(outputmodes), length(outputportnumbers),
            Nfrequencies),
        outputmode = outputmodes,
        outputport = outputportnumbers,
        freqindex=1:Nfrequencies,
    )
end

"""
    nodevariabletokeyed(nodevariable, outputmodes, nodenames)


# Examples
```jldoctest
julia> JosephsonCircuits.nodevariabletokeyed([1 2;3 4],[(0,),(1,)],["0","1","2"])
2-dimensional KeyedArray(NamedDimsArray(...)) with keys:
↓   outputmode ∈ 2-element Vector{Tuple{Int64}}
→   node ∈ 2-element Vector{String}
And data, 2×2 Matrix{Int64}:
          ("1")  ("2")
   (0,)    1      2
   (1,)    3      4
```
"""
function nodevariabletokeyed(nodevariable, outputmodes, nodenames)
    return  AxisKeys.KeyedArray(
        reshape(nodevariable, length(outputmodes), length(nodenames)-1),
        outputmode = outputmodes,
        node=nodenames[2:end])
end

"""
    nodevariabletokeyed(nodevariable, outputmodes, nodenames, inputmodes,
        inputportnumbers, w)

# Examples
```jldoctest
julia> JosephsonCircuits.nodevariabletokeyed([1 2;3 4;;;],[(0,),(1,)],["0","1"],[(0,),(1,)],[1],[1.0])
5-dimensional KeyedArray(NamedDimsArray(...)) with keys:
↓   outputmode ∈ 2-element Vector{Tuple{Int64}}
→   node ∈ 1-element Vector{String}
◪   inputmode ∈ 2-element Vector{Tuple{Int64}}
▨   inputport ∈ 1-element Vector{Int64}
▨   freqindex ∈ 1-element UnitRange{Int64}
And data, 2×1×2×1×1 Array{Int64, 5}:
[:, :, 1, 1, 1] ~ (:, :, (0,), 1, 1):
          ("1")
   (0,)    1
   (1,)    3

[:, :, 2, 1, 1] ~ (:, :, (1,), 1, 1):
          ("1")
   (0,)    2
   (1,)    4
```
"""
function nodevariabletokeyed(nodevariable, outputmodes, nodenames, inputmodes,
    inputportnumbers, w)

    return AxisKeys.KeyedArray(
        reshape(
            nodevariable,
            length(outputmodes),
            length(nodenames)-1,
            length(inputmodes),
            length(inputportnumbers),
            length(w),
        ),
        outputmode = outputmodes,
        node = nodenames[2:end],
        inputmode = inputmodes,
        inputport = inputportnumbers,
        freqindex=1:length(w),
    )
end

"""
    Snoisetokeyed(Snoise, inputmodes, components, outputmodes,
        outputportnumbers, w

# Examples
```jldoctest
julia> JosephsonCircuits.Snoisetokeyed([11 12;21 22;;;],[(0,)],["C1","C2"],[(0,)],[1,2],[1.0])
5-dimensional KeyedArray(NamedDimsArray(...)) with keys:
↓   inputmode ∈ 1-element Vector{Tuple{Int64}}
→   component ∈ 2-element Vector{String}
◪   outputmode ∈ 1-element Vector{Tuple{Int64}}
▨   outputport ∈ 2-element Vector{Int64}
▨   freqindex ∈ 1-element UnitRange{Int64}
And data, 1×2×1×2×1 Array{Int64, 5}:
[:, :, 1, 1, 1] ~ (:, :, (0,), 1, 1):
          ("C1")  ("C2")
   (0,)   11      21

[:, :, 1, 2, 1] ~ (:, :, (0,), 2, 1):
          ("C1")  ("C2")
   (0,)   12      22
```
"""
function Snoisetokeyed(Snoise, inputmodes, components, outputmodes,
    outputportnumbers, w)

    return AxisKeys.KeyedArray(
        reshape(
            Snoise,
            length(inputmodes),
            length(components),
            length(outputmodes),
            length(outputportnumbers),
            length(w),
        ),
        inputmode = inputmodes,
        component = components,
        outputmode = outputmodes,
        outputport = outputportnumbers,
        freqindex=1:length(w),
    )
end
