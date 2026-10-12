# The threads of the session, which the sweeps over frequencies, the
# batches of conditions, the network connections and the fits share: the
# chunks of work a batch is split into, run on tasks of their own with the
# first failure reported as the whole would meet it, and the number of
# BLAS threads, held at one while the threads factorize on their own.

# The conditions each task of a host batch steps. The conditions of a
# batch are independent of one another, so a host splits them across the
# threads of the session and steps the chunks at once: the whole step
# parallelizes that way, not only the assembly, the factorization and the
# solve, and the arrays of the batch are filled in place through views, so
# nothing is copied to put the batch back together. A device keeps the one
# chunk it has always had: its batch is already the parallelism, and
# driving one device from several tasks at once is not how the package
# uses it.
function batchchunks(backend, N::Integer)
    (backend isa CPU && N > 1) || return [1:N]
    nc = min(Base.Threads.nthreads(), N)
    nc > 1 || return [1:N]
    return [round(Int, (k - 1)*N/nc) + 1:round(Int, k*N/nc) for k in 1:nc]
end

# the error the chunks' tasks met, out of the failure of the tasks the
# threaded batch reports it as: the failed steps of the chunks as the
# first of them, with the conditions of every chunk which failed there,
# or else the first other error
function chunkerror(err::CompositeException)
    isempty(err.exceptions) && return err
    errs = [e for e in (chunkerror(e) for e in err.exceptions) if !(e isa ChunkAborted)]
    isempty(errs) && return err
    all(e -> e isa TransientStepError, errs) || return first(e for e in errs if !(e isa TransientStepError))
    step = minimum(e -> e.step, errs)
    at = [e for e in errs if e.step == step]
    return TransientStepError(step, first(at).time, sort!(reduce(vcat, [e.conditions for e in at])), first(at).cause)
end
chunkerror(err::TaskFailedException) = chunkerror(err.task.result)
chunkerror(err) = err

# The chunks of one `runchunks` share a flag, the step past which none
# goes on, kept in each chunk's task storage so that the stepping loops
# read it without an argument (`checkchunks`). A chunk whose step fails
# lowers it to that step: the others step on up to it, where a failure
# of theirs would still be the first, and stop there rather than at the
# end, so the error is the one the whole batch would meet. Any other
# error lowers it to zero, and a chunk stopped by it throws
# `ChunkAborted`, which the chunks' error leaves out.
struct ChunkAborted <: Exception end
const CHUNKFLAG = :transientchunkflag

# the stop of a chunk at a step past the flag its task holds, if it holds
# one; a step counts forward in a solve and a tangent and backward in an
# adjoint, whose errors all lower the flag to zero
function checkchunks(step::Integer)
    flag = get(task_local_storage(), CHUNKFLAG, nothing)
    flag isa Threads.Atomic{Int} && step > flag[] && throw(ChunkAborted())
    return nothing
end

# `f(k, chunk)` for every chunk of conditions, returned as a vector: on
# the calling task for one chunk, and otherwise each chunk on a task of
# its own across the threads of the session, which stop after a failure
# as the flag above says, a failure thrown as the chunks' error (see
# `chunkerror`). `f` is called once a chunk, and compiled for its own
# type; this is compiled once whatever closure it is given.
function runchunks(@nospecialize(f), chunks)
    length(chunks) == 1 && return Any[f(1, only(chunks))]
    each = Vector{Any}(undef, length(chunks))
    flag = Threads.Atomic{Int}(typemax(Int))
    try
        Base.Threads.@sync for (k, ch) in enumerate(chunks)
            Base.Threads.@spawn begin
                task_local_storage(CHUNKFLAG, flag)
                try
                    each[k] = f(k, ch)
                catch err
                    err isa ChunkAborted || Threads.atomic_min!(flag, err isa TransientStepError ? err.step : 0)
                    rethrow()
                end
            end
        end
    catch err
        throw(chunkerror(err))
    end
    return each
end

# The number of BLAS threads is a setting of the whole process, so the
# sweeps which run with one are counted: the first to start saves the
# setting and the last to finish restores it, which is right however many
# run at once.
const BLASLIMITED = Ref(0)
const BLASSAVED = Ref(1)
const BLASLIMITLOCK = ReentrantLock()

function withoneblasthread(f)
    lock(BLASLIMITLOCK) do
        if BLASLIMITED[] == 0
            BLASSAVED[] = BLAS.get_num_threads()
            BLAS.set_num_threads(1)
        end
        BLASLIMITED[] += 1
    end
    try
        return f()
    finally
        lock(BLASLIMITLOCK) do
            BLASLIMITED[] -= 1
            BLASLIMITED[] == 0 && BLAS.set_num_threads(BLASSAVED[])
        end
    end
end
