# Options shared by the benchmark scripts, included once per process.
#
# `--smoke` runs every case at a small size and with few samples, which is
# what the continuous integration job uses: it measures nothing, it only
# checks that each script still runs against the package's current
# interfaces. A label which is not an option names the revision in the
# output rows.
const BENCH_SMOKE = "--smoke" in ARGS
const BENCH_SAMPLES = BENCH_SMOKE ? 2 : 15
const BENCH_LABEL = let labels = filter(a -> !startswith(a, "--"), ARGS)
    isempty(labels) ? "current" : first(labels)
end

# the size of a case: the full size, or a small one under `--smoke`
benchsize(n::Int) = BENCH_SMOKE ? max(4, n ÷ 128) : n

# The timer of every script: `f(x)` timed on its first call, called once
# more to warm up, and timed `BENCH_SAMPLES` times. Returns the first
# call's `@timed` record, its compile time, which `@timed` reports from
# Julia 1.11 on and is NaN before, and the samples.
function benchtime(f, x)
    @nospecialize f x
    GC.gc()
    first = @timed Base.invokelatest(f, x)
    Base.invokelatest(f, x)
    samples = [@timed Base.invokelatest(f, x) for _ in 1:BENCH_SAMPLES]
    compile = hasproperty(first, :compile_time) ? first.compile_time : NaN
    return first, compile, samples
end
