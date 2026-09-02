"""
    DocstringCheck

Compare the signature lines at the top of each docstring in a set of Julia
source files with the definition the docstring is attached to, and report
where the two disagree.

The check is textual: a docstring is a line holding only `\"\"\"`, its
signature block is the run of four space indented lines at its top, and its
definition is the first non blank, non comment line after the closing
`\"\"\"`. A signature or definition spanning several lines is joined until
its parentheses balance; a multi line string literal opened on a code line
is skipped. Positional argument names are compared in order, keyword
argument names as sets. Argument names are the text before any `=` default
and `::` type. A keyword splat (`kwargs...`) on either side stands for
keywords documented elsewhere: with one in the definition, documented
keywords it forwards are accepted, and with one in the docstring,
definition keywords the docstring omits are accepted.

Reported issues are `Issue`s of kind `:name` (no signature line names the
definition), `:args` (a signature names it but its arguments differ) or
`:indent` (the signature block is indented by other than four spaces, so
Documenter renders it as prose). A docstring with no signature block is not
checked.

Use [`checkdocstrings`](@ref) on a package root, or run this file as a
script with the root as its argument.
"""
module DocstringCheck

export Issue, checkdocstrings, checkfile, format

"""
    Issue(file, line, kind, name, detail)

One disagreement between a docstring and its definition. `kind` is `:name`,
`:args` or `:indent`; `line` is the last line of the docstring's signature
block.
"""
struct Issue
    file::String
    line::Int
    kind::Symbol
    name::String
    detail::String
end

# The definition keywords whose docstrings are not checked: they carry no
# argument list (or, for `const`, an arbitrary right hand side).
const SKIPPED_HEADERS = ("struct ", "mutable struct ", "abstract type ",
    "const ", "module ", "macro ")

# Prefixes stripped from a definition before its name and arguments are
# read. They are stripped repeatedly, so `@inline function f(x)` reads as
# `f(x)`.
const HEADER_PREFIXES = ("@inline ", "@noinline ", "@kernel ", "function ")

# --- parsing ---------------------------------------------------------------

"""
    argnamelist(s)

The bare argument names of the comma separated list `s`: split at top level
commas, then each part is cut at its first `=` and its first `::`, stripped,
and stripped of `...`. An unnamed argument (`::Type{T}`) gives `""`.
"""
function argnamelist(s::AbstractString)
    parts = String[]
    depth = 0
    cur = IOBuffer()
    for c in s
        if c in ('(', '[', '{')
            depth += 1
        elseif c in (')', ']', '}')
            depth -= 1
        end
        if c == ',' && depth == 0
            push!(parts, String(take!(cur)))
        else
            write(cur, c)
        end
    end
    push!(parts, String(take!(cur)))
    names = String[]
    for a in parts
        a = strip(a)
        isempty(a) && continue
        a = first(split(a, '='; limit = 2))
        a = first(split(a, "::"; limit = 2))
        push!(names, String(strip(a)))
    end
    return names
end

# split a name list into the plain names and whether it held a splat
# (`kwargs...`), which stands for keywords documented elsewhere
function splitsplat(names::Vector{String})
    plain = String[]
    splat = false
    for a in names
        if endswith(a, "...")
            splat = true
        else
            push!(plain, a)
        end
    end
    return plain, splat
end

"""
    argnames(header)

The positional and keyword argument names of the call or definition
`header`, read from its first balanced parenthesized list, or `nothing` when
it has none. The list is split at its first `;` into the two groups.
"""
function argnames(header::AbstractString)
    start = findfirst('(', header)
    start === nothing && return nothing
    depth = 0
    stop = nothing
    for k in start:lastindex(header)
        isvalid(header, k) || continue
        c = header[k]
        if c == '('
            depth += 1
        elseif c == ')'
            depth -= 1
            if depth == 0
                stop = k
                break
            end
        end
    end
    stop === nothing && return nothing
    inner = header[nextind(header, start):prevind(header, stop)]
    semi = findfirst(';', inner)
    pos, kw = semi === nothing ? (inner, "") :
        (inner[1:prevind(inner, semi)], inner[nextind(inner, semi):end])
    return argnamelist(pos), argnamelist(kw)
end

# the name a definition or signature line starts with, without a module
# prefix or a macro sigil, or `nothing`
function leadingname(s::AbstractString)
    m = match(r"^[\w.!@]+", s)
    m === nothing && return nothing
    name = last(split(m.match, '.'))
    return lstrip(name, '@')
end

stripwhere(s::AbstractString) = String(strip(replace(s, r"\bwhere\b.*$" => "")))

function stripprefixes(s::AbstractString)
    s = String(strip(s))
    changed = true
    while changed
        changed = false
        for p in HEADER_PREFIXES
            if startswith(s, p)
                s = String(strip(s[length(p)+1:end]))
                changed = true
            end
        end
    end
    return s
end

# positional lists agree when they have the same length and each pair of
# names is equal, an unnamed argument on either side matching any name
function positionalsagree(doc::Vector{String}, code::Vector{String};
        wildcardunnamed::Bool = true)
    length(doc) == length(code) || return false
    for (a, b) in zip(doc, code)
        a == b && continue
        wildcardunnamed && (isempty(a) || isempty(b)) && continue
        return false
    end
    return true
end

# a Python style list repr, so that the output can be diffed against the
# original Python checker's
pyrepr(v::Vector{String}) = "[" * join(("'" * x * "'" for x in v), ", ") * "]"

# --- the check -------------------------------------------------------------

"""
    checkfile(path; wildcardunnamed = true) -> Vector{Issue}

Check every docstring in the file at `path`. With `wildcardunnamed` an
unnamed positional argument in either the docstring or the definition
matches any name; without it the two must be spelled identically.
"""
function checkfile(path::AbstractString; wildcardunnamed::Bool = true)
    lines = readlines(path)
    n = length(lines)
    issues = Issue[]
    i = 1
    while i <= n
        # a code line which opens a multi line string literal (an odd number
        # of `"""` on it, as in `const text = """`) is skipped to its close,
        # so the closing `"""` is not mistaken for a docstring opener
        if strip(lines[i]) != "\"\"\"" && isodd(count("\"\"\"", lines[i]))
            i += 1
            while i <= n && !occursin("\"\"\"", lines[i])
                i += 1
            end
            i += 1
            continue
        end
        # a docstring opens on a line holding only `"""`, unless the
        # previous line ends in `=`, which makes it a string constant
        if !(strip(lines[i]) == "\"\"\"" &&
                (i == 1 || !endswith(strip(lines[i-1]), '=')))
            i += 1
            continue
        end
        j = i + 1
        while j <= n && isempty(strip(lines[j]))
            j += 1
        end
        # the signature block: consecutive indented lines, joined while their
        # parentheses are unbalanced. Documenter needs four spaces for a code
        # block, so a shallower indent is reported.
        sigs = String[]
        cur = ""
        indentok = true
        while j <= n && startswith(lines[j], " ")
            startswith(lines[j], "    ") || (indentok = false)
            l = String(strip(lines[j]))
            if !isempty(cur) && count('(', cur) == count(')', cur)
                push!(sigs, cur)
                cur = l
            else
                cur = String(strip(cur * " " * l))
            end
            j += 1
        end
        isempty(cur) || push!(sigs, cur)
        reportline = j - 1   # the last line of the signature block
        if !indentok
            push!(issues, Issue(path, reportline, :indent,
                isempty(sigs) ? "" : first(sigs),
                "signature block indented by fewer than four spaces"))
        end
        # the closing `"""`, then the first non blank, non comment line
        k = j
        while k <= n && strip(lines[k]) != "\"\"\""
            k += 1
        end
        m = k + 1
        while m <= n && (isempty(strip(lines[m])) || startswith(strip(lines[m]), '#'))
            m += 1
        end
        if m > n || isempty(sigs)
            i = k + 1
            continue
        end
        hdr = String(strip(lines[m]))
        q = m
        while count('(', hdr) > count(')', hdr) && q + 1 <= n
            q += 1
            hdr *= " " * strip(lines[q])
        end
        if any(p -> startswith(hdr, p), SKIPPED_HEADERS)
            i = k + 1
            continue
        end
        hdr2 = stripwhere(stripprefixes(hdr))
        dargs = argnames(hdr2)
        dname = leadingname(hdr2)
        if dargs === nothing || dname === nothing
            i = k + 1
            continue
        end
        dpos, dkwall = dargs
        dpos = [replace(x, "..." => "") for x in dpos]
        dkw, dsplat = splitsplat(dkwall)
        ok = false
        best = nothing
        for s in sigs
            s2 = stripwhere(s)
            sname = leadingname(s2)
            sname === nothing && continue
            sargs = argnames(s2)
            sargs === nothing && continue
            sname == dname || continue
            spos = [replace(x, "..." => "") for x in sargs[1]]
            skw, ssplat = splitsplat(sargs[2])
            # a splat on either side stands for keywords documented elsewhere:
            # the docstring may then omit definition keywords, and the
            # definition may forward documented ones
            missing = dsplat ? String[] : [x for x in skw if !(x in dkw)]
            extra = ssplat ? String[] : [x for x in dkw if !(x in skw)]
            best = (s, spos, missing, extra)
            if positionalsagree(spos, dpos; wildcardunnamed) &&
                    isempty(missing) && isempty(extra)
                ok = true
                break
            end
        end
        if !ok
            if best === nothing
                push!(issues, Issue(path, reportline, :name,
                    isempty(sigs) ? "" : first(sigs), hdr))
            else
                s, spos, missing, extra = best
                parts = String[]
                positionalsagree(spos, dpos; wildcardunnamed) ||
                    push!(parts, "pos doc=" * pyrepr(spos) * " code=" * pyrepr(dpos))
                isempty(missing) || push!(parts, "doc-only kw=" * join(missing, ","))
                isempty(extra) || push!(parts, "undocumented kw=" * join(extra, ","))
                isempty(parts) || push!(issues, Issue(path, reportline, :args,
                    dname, join(parts, "; ")))
            end
        end
        i = k + 1
    end
    return issues
end

"""
    checkdocstrings(root; dirs = ("src", "ext"), wildcardunnamed = true,
        ignore = _ -> false) -> Vector{Issue}

Check every `.jl` file under the given subdirectories of `root`. `ignore`
is a predicate on `Issue`s which are known and accepted; those are dropped.
"""
function checkdocstrings(root::AbstractString; dirs = ("src", "ext"),
        wildcardunnamed::Bool = true, ignore = _ -> false)
    issues = Issue[]
    for d in dirs
        dir = joinpath(root, d)
        isdir(dir) || continue
        for f in sort(readdir(dir; join = true))
            endswith(f, ".jl") || continue
            for iss in checkfile(f; wildcardunnamed)
                ignore(iss) || push!(issues, iss)
            end
        end
    end
    return issues
end

"""
    format(issue::Issue, root = "")

One line describing `issue`, in the form `file:line | KIND | name | detail`.
"""
function format(iss::Issue, root::AbstractString = "")
    f = isempty(root) ? iss.file : relpath(iss.file, root)
    kind = iss.kind === :name ? "NAME" : iss.kind === :indent ? "INDENT" : "ARGS"
    return "$(f):$(iss.line) | $(kind) | $(iss.name) | $(iss.detail)"
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    root = isempty(ARGS) ? pwd() : ARGS[1]
    wildcard = !("--strict" in ARGS)
    issues = DocstringCheck.checkdocstrings(root; wildcardunnamed = wildcard)
    for iss in issues
        println(DocstringCheck.format(iss, root))
    end
    println(length(issues), " issues")
end
