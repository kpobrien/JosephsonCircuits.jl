
# Deprecated `connectS` and `connectS!` intra and interconnection functions.
# These are deprecated in order to allow the `intraconnectS` and
# `interconnectS` functions to accept noise correlation matrices without the
# ambiguity of whether the second matrix is a second scattering
# parameter matrix or a noise correlation matrix.
function connectS(Sa::AbstractArray{T,N}, k::Int, l::Int;
    nbatches::Int = Base.Threads.nthreads()) where {T,N}
    Base.depwarn("connectS(Sa::AbstractArray, k::Int, l::Int)` is deprecated, use `intraconnectS(Sa, k, l)` instead.", :connectS; force=true)
    return intraconnectS(Sa, k, l; nbatches = nbatches)
end

function connectS(Sa::AbstractArray{T,N}, Sb::AbstractArray{T,N}, k::Int, l::Int;
    nbatches::Int = Base.Threads.nthreads()) where {T,N}
    Base.depwarn("connectS(Sa::AbstractArray, Sb::AbstractArray, k::Int, l::Int)` is deprecated, use `interconnectS(Sa, Sb, k, l)` instead.", :connectS; force=true)
    return interconnectS(Sa, Sb, k, l; nbatches = nbatches)
end

function connectS!(Sout, Sa, k::Int, l::Int;
    nbatches::Int = Base.Threads.nthreads())
    Base.depwarn("connectS!(Sout, Sa, k::Int, l::Int)` is deprecated, use `intraconnectS!(Sout, Sa, k, l)` instead.", :connectS!; force=true)
    return intraconnectS!(Sout, Sa, k, l; nbatches = nbatches)
end

function connectS!(Sout, Sa, Sb, k::Int, l::Int;
    nbatches::Int = Base.Threads.nthreads())
    Base.depwarn("connectS!(Sout, Sa, Sb, k::Int, l::Int)` is deprecated, use `interconnectS!(Sout, Sa, Sb, k, l)` instead.", :connectS!; force=true)
    return interconnectS!(Sout, Sa, Sb, k, l; nbatches = nbatches)
end