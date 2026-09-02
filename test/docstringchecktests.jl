# Tests of the docstring signature checker in docstringcheck.jl: synthetic
# files with known agreements and disagreements, then the package itself.

using Test
include("docstringcheck.jl")
using .DocstringCheck

# write `src` as a source file in a temporary package root and check it
function checksource(src::AbstractString; kwargs...)
    mktempdir() do root
        mkdir(joinpath(root, "src"))
        write(joinpath(root, "src", "t.jl"), src)
        return checkdocstrings(root; kwargs...)
    end
end

@testset "DocstringCheck" begin

    @testset "argument name parsing" begin
        @test DocstringCheck.argnamelist("a, b::Int, c = 1, d::T = 2, e...") ==
            ["a", "b", "c", "d", "e..."]
        @test DocstringCheck.splitsplat(["c", "kwargs..."]) == (["c"], true)
        @test DocstringCheck.splitsplat(["c"]) == (["c"], false)
        @test DocstringCheck.argnamelist("::Type{T}, n") == ["", "n"]
        @test DocstringCheck.argnamelist("f(x, y), g = (1, 2), h") ==
            ["f(x, y)", "g", "h"]
        @test DocstringCheck.argnamelist("") == String[]
        @test DocstringCheck.argnames("f(a, b; c = 1, kwargs...)") ==
            (["a", "b"], ["c", "kwargs..."])
        @test DocstringCheck.argnames("f(a)") == (["a"], String[])
        @test DocstringCheck.argnames("nothing here") === nothing
        @test DocstringCheck.argnames("f(x::NTuple{N,Int}, y) where N") ==
            (["x", "y"], String[])
        @test DocstringCheck.leadingname("Base.show(io, x)") == "show"
        @test DocstringCheck.leadingname("@params a b") == "params"
        @test DocstringCheck.stripprefixes("@inline function f(x)") == "f(x)"
        @test DocstringCheck.stripprefixes("@kernel function k!(a)") == "k!(a)"
        @test DocstringCheck.stripprefixes("function  g(x)") == "g(x)"
        @test DocstringCheck.stripwhere("f(x::T) where {T<:Real}") == "f(x::T)"
    end

    @testset "agreeing docstrings report nothing" begin
        @test isempty(checksource("""
            \"\"\"
                f(a, b; c = 1)

            Doc.
            \"\"\"
            function f(a, b; c = 1)
                a
            end

            \"\"\"
                g(x::Int; y = 2, kwargs...)
                g(x::Float64; y = 2)

            Two methods.
            \"\"\"
            function g(x::Float64; y = 2)
                x
            end

            \"\"\"
                h(a,
                    b)

            A signature spanning two lines.
            \"\"\"
            function h(a,
                    b)
                a
            end

            \"\"\"
                Base.show(io::IO, p::Port)

            A qualified name.
            \"\"\"
            function Base.show(io::IO, p::Port)
            end

            \"\"\"
                k(x::T) where T

            A where clause on both sides.
            \"\"\"
            @inline function k(x::T) where {T<:Real}
                x
            end
            """))
    end

    @testset "keyword mismatches" begin
        issues = checksource("""
            \"\"\"
                f(a; b = 1)

            Doc.
            \"\"\"
            function f(a; c = 1)
                a
            end
            """)
        @test length(issues) == 1
        @test issues[1].kind === :args
        @test issues[1].name == "f"
        @test issues[1].line == 2
        @test occursin("doc-only kw=b", issues[1].detail)
        @test occursin("undocumented kw=c", issues[1].detail)
        # a splat in the definition may forward documented keywords, and a
        # splat in the docstring may stand for undocumented definition ones
        @test isempty(checksource("""
            \"\"\"
                f(a; b = 1, c = 2)
            \"\"\"
            f(a; b = 1, kw...) = a
            """))
        @test isempty(checksource("""
            \"\"\"
                g(a; b = 1, kwargs...)
            \"\"\"
            g(a; b = 1, c = 2, d = 3) = a
            """))
        # but a splat does not excuse a keyword the definition lacks when
        # the definition has none
        @test length(checksource("""
            \"\"\"
                h(a; b = 1, kwargs...)
            \"\"\"
            h(a; c = 1) = a
            """)) == 1
    end

    @testset "positional mismatches" begin
        issues = checksource("""
            \"\"\"
                f(a, b)

            Doc.
            \"\"\"
            function f(a, b, c)
            end
            """)
        @test length(issues) == 1 && issues[1].kind === :args
        @test occursin("pos doc=['a', 'b'] code=['a', 'b', 'c']", issues[1].detail)
        # an unnamed argument matches any documented name by default, and
        # nothing in strict mode
        src = """
            \"\"\"
                f(m::AbstractDCModel, n)
            \"\"\"
            f(::OpenDC, n) = n
            """
        @test isempty(checksource(src))
        @test length(checksource(src; wildcardunnamed = false)) == 1
    end

    @testset "name mismatches" begin
        issues = checksource("""
            \"\"\"
                oldname(a)

            Doc.
            \"\"\"
            function newname(a)
            end
            """)
        @test length(issues) == 1
        @test issues[1].kind === :name
        @test issues[1].name == "oldname(a)"
        @test startswith(issues[1].detail, "function newname(a)")
    end

    @testset "skipped definitions" begin
        # structs, constants, modules, macros, and string constants are not
        # checked; neither is a docstring with no signature block
        @test isempty(checksource("""
            \"\"\"
                Foo(a, b)

            A struct.
            \"\"\"
            struct Foo
                a
                b
            end

            \"\"\"
                bar

            A constant.
            \"\"\"
            const bar = 1

            const text = \"\"\"
                notasignature(x)
            \"\"\"
            function unrelated(y)
            end

            \"\"\"
            No signature block at all.
            \"\"\"
            function alsofine(z)
            end

            \"\"\"
                trailing(a)

            A docstring at the end of the file.
            \"\"\"
            """))
    end

    @testset "misindented signature block" begin
        issues = checksource("""
            \"\"\"
              f(a)

            Two spaces: Documenter shows this as prose.
            \"\"\"
            function f(a)
            end
            """)
        @test length(issues) == 1
        @test issues[1].kind === :indent
        @test issues[1].line == 2
        @test occursin("INDENT", format(issues[1]))
    end

    @testset "comments between docstring and definition" begin
        issues = checksource("""
            \"\"\"
                f(a)
            \"\"\"
            # a comment, and

            # a blank line, are skipped
            function f(b)
            end
            """)
        @test length(issues) == 1 && issues[1].kind === :args
    end

    @testset "ignore predicate and formatting" begin
        src = """
            \"\"\"
                f(a; b = 1)
            \"\"\"
            f(a; c = 1) = a
            """
        @test isempty(checksource(src; ignore = iss -> iss.name == "f"))
        issues = checksource(src)
        line = format(issues[1])
        @test occursin("t.jl:2 | ARGS | f | ", line)
    end

    @testset "the package itself" begin
        root = joinpath(@__DIR__, "..")
        # accepted: docstrings which name an argument differently from the
        # definition for readability, kernels documented by name only, and
        # the deprecated or removed keywords the solver entry points still
        # accept but document in prose rather than in their signature lines
        deprecated = ("switchofflinesearchtol", "alphamin", "returnZ",
            "returnZadjoint", "returnZsensitivity", "returnZsensitivityadjoint",
            "sensitivitypairs", "sensitivityblockpairs",
            "nsensitivityparameters", "sensitivitylabels")
        function accepted(iss)
            occursin("pos doc=", iss.detail) && return true
            iss.kind === :name && return true
            m = match(r"undocumented kw=([\w,]+)", iss.detail)
            m !== nothing && all(in(deprecated), split(m.captures[1], ',')) &&
                return true
            return false
        end
        issues = checkdocstrings(root; ignore = accepted)
        for iss in issues
            println(format(iss, root))
        end
        @test isempty(issues)
    end
end
