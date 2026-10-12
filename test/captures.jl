using JosephsonCircuits
using Test

# The variables which closures capture boxed in the methods the package
# and its loaded extensions define, methods of other modules' functions
# included, from their lowered code, and the number of methods read. A
# variable which a closure captures and which is assigned more than once,
# or a variable of the enclosing function which a closure assigns, is held
# in a box: each use of it is a dynamic call on a value of unknown type,
# and a closure meant to have a local of its own writes the enclosing
# function's.
function boxedcaptures()
    owners = Set{Module}(m for m in Base.loaded_modules_array()
        if startswith(string(nameof(m)), "JosephsonCircuits"))
    seen = Base.IdSet{Method}()
    boxed = String[]
    for mod in Base.loaded_modules_array(), n in names(mod; all = true)
        isdefined(mod, n) || continue
        f = getglobal(mod, n)
        (f isa Function || f isa Type) || continue
        ms = try
            methods(f)
        catch
            continue
        end
        for m in ms
            (m.module in owners && !(m in seen)) || continue
            push!(seen, m)
            ci = Base.uncompressed_ast(m)
            vars = [ci.slotnames[st.args[1].id] for st in ci.code if st isa Expr &&
                st.head === :(=) && st.args[1] isa Core.SlotNumber && occursin("Core.Box", string(st.args[2]))]
            isempty(vars) || push!(boxed, "$(m.file):$(m.line) $(m.name): $(join(vars, ", "))")
        end
    end
    return boxed, length(seen)
end

@testset "closures capture no variable boxed" begin
    boxed, scanned = boxedcaptures()
    # the scan read the package's methods
    @test scanned > 1000
    @test boxed == String[]
end
