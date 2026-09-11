"""
    JosephsonCircuitsXicToolsExt

The `wrspice` command of XicTools_jll as the default executable of
[`JosephsonCircuits.wrspice_cmd`](@ref): loading XicTools_jll lets
`WRspice()` and `spice_run` find WRSPICE without a path, on the
platforms the artifact supports.
"""
module JosephsonCircuitsXicToolsExt

import JosephsonCircuits
import XicTools_jll

function __init__()
    if XicTools_jll.is_available()
        JosephsonCircuits.wrspicedefaultcmd[] = XicTools_jll.wrspice()
    end
    return nothing
end

end
