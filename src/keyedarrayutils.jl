

"""
    Stokeyed(S, outputmodes, outputportnumbers, inputmodes, inputportnumbers, w)

Convert a scattering parameter array to a keyed array. 

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

function nodevariabletokeyed(nodevariable, outputmodes, nodes)
    return  AxisKeys.KeyedArray(
        reshape(nodevariable, length(outputmodes), length(nodes)),
        outputmode = outputmodes,
        node=nodes)
end


function nodevariabletokeyed(nodevariable, outputmodes, nodes, inputmodes,
    inputportnumbers, w)

    return AxisKeys.KeyedArray(
        reshape(
            nodevariable,
            length(outputmodes),
            length(nodes),
            length(inputmodes),
            length(inputportnumbers),
            length(w),
        ),
        outputmode = outputmodes,
        node = nodes,
        inputmode = inputmodes,
        inputport = inputportnumbers,
        freqindex=1:length(w),
    )
end

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
