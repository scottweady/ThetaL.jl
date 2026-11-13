# ThetaL's code map: which functions of src/ each feature runs. Every setting (a scheme, a kind of velocity or symbol,
# an integrator feature, a driver, a curve utility) is one small run, changed from the base case (curve shortening
# of a star with ETDRK2, two steps of `advance!`) in that setting only, run under line coverage; the lines it
# executes are mapped back to the top-level function definitions of src/. Written as data for
# src/assets/codemap/codemap.js, which draws the docs page "Code map".
#   julia -t 8 --project=docs docs/codemap.jl [covdir]          (from the package root; a few minutes)
# The script runs itself again in a child Julia with coverage on (no package image, so the child compiles ThetaL).
# The child's tracefiles (<setting>.info) go to `covdir` if given, and a `covdir` that already holds them skips the run.
const ROOT = dirname(@__DIR__)                       # the package root
const SRC = joinpath(ROOT, "src")
const OUT = joinpath(@__DIR__, "src", "assets", "codemap", "codemap_data.js")

# ---- the settings: (id, label, group, what it is, the run) ----------------------------------------------------
# `parent` is the scheme an option is run with (the page compares it with that scheme); "" for the schemes. `schemes`:
# the schemes the option combines with (all of them, here).
const BOTH = ["etdrk2", "expeuler"]
const SETTINGS = [
    (id="etdrk2", label="ETDRK2", group="Scheme", base=true, parent="", schemes=String[],
     about="ETDRK2 (the default): curve shortening of a star, init and two steps of advance! (the base case)",
     run=:(begin integ = init(CurveProblem([star], mcf; symbol=q -> -q^2), ETDRK2()); advance!(integ, 1e-3); advance!(integ, 2e-3) end)),
    (id="expeuler", label="ExpEuler", group="Scheme", base=false, parent="", schemes=String[],
     about="ExpEuler: exponential Euler, one evaluation per step",
     run=:(begin integ = init(CurveProblem([star], mcf; symbol=q -> -q^2), ExpEuler()); advance!(integ, 1e-3); advance!(integ, 2e-3) end)),
    (id="nosymbol", label="No symbol", group="Velocity", base=false, parent="etdrk2", schemes=BOTH,
     about="symbol = nothing: no stiff part (the explicit midpoint rule)",
     run=:(begin integ = init(CurveProblem([star], (Xs, ev; warm) -> [fill(0.1, length(X.φ)) for X in Xs]), ETDRK2()); advance!(integ, 1e-3) end)),
    (id="statesymbol", label="State-dependent symbol", group="Velocity", base=false, parent="etdrk2", schemes=BOTH,
     about="symbol(q, X) read from the curve X (here the largest curvature), frozen per step",
     run=:(begin integ = init(CurveProblem([star], mcf; symbol=(q, X) -> -maximum(abs, curvature(X)) * q^2), ETDRK2()); advance!(integ, 1e-3) end)),
    (id="flow", label="Flow field", group="Velocity", base=false, parent="etdrk2", schemes=BOTH,
     about="U = u·n of a velocity field u (a rotation about an outside point): velocity_components, the rigid part exact",
     run=:(begin integ = init(CurveProblem([star], (Xs, ev; warm) -> [velocity_components(X, [(-x[2] + 2, x[1]) for x in curve_points(X)])[1] for X in Xs])); advance!(integ, 0.1) end)),
    (id="norigid", label="No rigid split", group="Velocity", base=false, parent="etdrk2", schemes=BOTH,
     about="rigid = :none: the plain Hou–Lowengrub–Shelley step",
     run=:(begin integ = init(CurveProblem([star], mcf; symbol=q -> -q^2, rigid=:none), ETDRK2()); advance!(integ, 1e-3) end)),
    (id="coupled", label="Several curves", group="Velocity", base=false, parent="etdrk2", schemes=BOTH,
     about="two curves coupled through a mean field, with tags and a context (a named-tuple result)",
     run=:(begin integ = init(CurveProblem([star, θ -> (4 + 0.5cos(θ), 0.5sin(θ))], (Xs, ev; warm) -> (; U=[2π * length(Xs) / sum(perimeter, Xs) .- curvature(X) for X in Xs], ctx=ev.t); symbol=q -> -q^2, tags=[:a, :b]), ETDRK2()); advance!(integ, 1e-3) end)),
    (id="edits", label="Curve edits", group="Integrator", base=false, parent="etdrk2", schemes=BOTH,
     about="add_curve!, remove_curve!, set_curve!, offset_curve!, double_grid!, halve_grid!, curve_index",
     run=:(begin integ = init(CurveProblem([star], mcf; symbol=q -> -q^2)); advance!(integ, 1e-3)
                 id = add_curve!(integ, θ -> (4 + 0.5cos(θ), 0.5sin(θ))); double_grid!(integ, 1); halve_grid!(integ, 1)
                 offset_curve!(integ, curve_index(integ, id), 0.1); set_curve!(integ, 2, θ -> (4 + 0.4cos(θ), 0.4sin(θ)))
                 advance!(integ, 2e-3); remove_curve!(integ, 2); advance!(integ, 3e-3) end)),
    (id="reject", label="Rejection", group="Integrator", base=false, parent="etdrk2", schemes=BOTH,
     about="a StepFailure from the velocity, the error estimate against errtol, validate, and the retries",
     run=:(begin integ = init(CurveProblem([star], (Xs, ev; warm) -> (ev.kind === :stage && ev.attempt == 1 && throw(StepFailure(:test)); mcf(Xs, ev; warm)); symbol=q -> -q^2, validate=Xs -> true))
                 advance!(integ, 1e-3); advance!(integ, 1e-3); advance!(integ, 1.0; errtol=1e-12); step!(integ, 1e-4) end)),
    (id="checkpoint", label="Checkpoint", group="Integrator", base=false, parent="etdrk2", schemes=BOTH,
     about="checkpoint(integ) as plain data, and init(problem, scheme; from) to restart",
     run=:(begin prob = CurveProblem([star], mcf; symbol=q -> -q^2); integ = init(prob); advance!(integ, 1e-3)
                 again = init(prob, ETDRK2(); from=checkpoint(integ)); advance!(again, 2e-3) end)),
    (id="autofilter", label="Automatic filter", group="Integrator", base=false, parent="etdrk2", schemes=BOTH,
     about="filter_cutoff = :auto: the cutoff set from noise_cutoff after the first evaluation",
     run=:(begin integ = init(CurveProblem([star], mcf; symbol=q -> -q^2), ETDRK2(filter_cutoff=:auto)); advance!(integ, 1e-3) end)),
    (id="evolve", label="evolve", group="Driver", base=false, parent="etdrk2", schemes=BOTH,
     about="evolve(curve, U, tspan): one call, saved times, sol(t), save_curves",
     run=:(begin sol = evolve(star, (X, t) -> -curvature(X), (0, 2e-3); dt=1e-3, symbol=q -> -q^2, saveat=1e-3); sol(1e-3); save_curves(IOBuffer(), sol) end)),
    (id="solve", label="solve", group="Driver", base=false, parent="etdrk2", schemes=BOTH,
     about="solve(problem, scheme; tspan, dt, saveat, callback) on a CurveProblem",
     run=:(solve(CurveProblem([star], mcf; symbol=q -> -q^2), ETDRK2(); tspan=(0, 2e-3), dt=(Xs, t) -> 1e-3, callback=integ -> nothing))),
    (id="resolution", label="Resolution", group="Driver", base=false, parent="etdrk2", schemes=BOTH,
     about="N doubled past tol (a star on 64 points), and halved with coarsen = true (a circle on 512 points)",
     run=:(begin evolve(from_curve(star, 64), (X, t) -> zeros(length(X.φ)), (0, 0.2); dt=0.1)
                 evolve(ThetaLState(zeros(512), 2π, ThetaL.Vec2(0, -1)), (X, t) -> fill(0.1, length(X.φ)), (0, 0.12); dt=0.01, coarsen=true) end)),
    (id="construct", label="Curves from data", group="Curves", base=false, parent="etdrk2", schemes=BOTH,
     about="from_curve (fixed or automatic N), from_arclength_samples, curve_function",
     run=:(begin X = from_curve(star, 128); from_curve(θ -> (cos(θ), 2sin(θ)); tol=1e-8); from_arclength_samples(curve_points(X)); curve_function(X)(0.3) end)),
    (id="measures", label="Measures", group="Curves", base=false, parent="etdrk2", schemes=BOTH,
     about="min_distance, curve_integral, offset_curve, remove_area_rate, frame_velocity, rigid_part, geometry",
     run=:(begin X = from_curve(star, 128); U = -curvature(X); t = from_curve(θ -> (4 + cos(θ), sin(θ)), 64)
                 min_distance([X, t]); curve_integral(X, U); offset_curve(X, 0.1); remove_area_rate([X, t], [U, zeros(64)]; per_curve=false)
                 frame_velocity(X, U); rigid_part(X, U; mode=:translation); area(X); centroid(X); curve_moments(X)
                 arclength_derivative(X, U); velocity_vectors(X, U, U); spectral_tail(X) end)),
]

# ---- child: run the settings under coverage, one LCOV tracefile each -------------------------------------------
if Base.JLOptions().code_coverage != 0
    outdir = ARGS[1]
    using ThetaL
    Core.eval(Main, :(star = θ -> (1 + 0.3cos(5θ)) .* (cos(θ), sin(θ))))
    Core.eval(Main, :(mcf(Xs, ev; warm) = [-curvature(X) for X in Xs]))
    for X in SETTINGS
        sec = @elapsed begin
            ccall(:jl_clear_coverage_data, Cvoid, ())
            Core.eval(Main, X.run)
            ccall(:jl_write_coverage_data, Cvoid, (Cstring,), joinpath(outdir, X.id * ".info"))
        end
        println(rpad(X.id, 12), round(sec; digits=1), " s"); flush(stdout)
    end
    exit()
end

# ---- parent: the top-level functions of src/ and their line ranges --------------------------------------------
# A function is a name in a file (all its methods there); its lines run from its definition to the next top-level
# statement, so a docstring belongs to the function it documents.
_unwrap(ex) = ex isa Expr && ex.head === :macrocall ? _unwrap(ex.args[end]) : ex
_docstring(ex) = ex isa Expr && ex.head === :macrocall && ex.args[1] == GlobalRef(Core, Symbol("@doc")) &&
                 ex.args[3] isa String ? ex.args[3] : nothing
const GENERIC = Set([:size, :eltype, :length, :show, :getindex, :iterate, :copy, :convert, :(==)])
const BASETYPES = Set([:AbstractVector, :AbstractMatrix, :AbstractArray, :Vector, :Matrix, :Number, :Real, :Integer,
                       :Int, :Float64, :IO, :MIME, :Function, :Any, :Tuple, :NamedTuple])
_tname(t) = t isa Symbol ? t : t isa Expr && t.head === :curly ? _tname(t.args[1]) :
            t isa Expr && t.head === :. ? t.args[end].value : nothing
function _argtype(sig)
    for a in sig.args[2:end]
        a isa Expr && a.head === :parameters && continue
        a isa Expr && a.head === :kw && (a = a.args[1])
        a isa Expr && a.head === :(::) || continue
        t = _tname(a.args[end])
        t === nothing || t in BASETYPES || return t
    end
    nothing
end
function _fname(ex)
    ex = _unwrap(ex)
    ex isa Expr && ex.head in (:function, :(=)) || return nothing
    sig = ex.args[1]
    while sig isa Expr && sig.head in (:where, :(::)) && !(sig.head === :(::) && length(sig.args) == 1)
        sig = sig.args[1]
    end
    sig isa Expr && sig.head === :call || return nothing
    f = sig.args[1]
    f isa Expr && f.head === :curly && (f = f.args[1])
    f isa Expr && f.head === :. && (f = f.args[end].value)                          # CommonSolve.init, Base.show
    f isa Expr && f.head === :(::) && return "(::" * string(_tname(f.args[end])) * ")"   # (op::T)(x)
    f isa Symbol || return nothing
    t = f in GENERIC ? _argtype(sig) : nothing
    t === nothing ? string(f) : string(f, "(::", t, ")")
end
function _summary(doc)
    doc === nothing && return ""
    lines = split(strip(doc), '\n')
    k = findfirst(l -> !startswith(l, "    ") && !isempty(strip(l)), lines)       # skip the signature block
    k === nothing && return ""
    para = String[]
    for l in lines[k:end]
        isempty(strip(l)) && break
        push!(para, strip(l))
    end
    X = join(para, " ")
    m = match(r"^(.*?\.)(\X|$)", X)
    m === nothing ? X : m.captures[1]
end
function functions_of(path)
    text = read(path, String)
    nlines = count(==('\n'), text) + 1
    top = Meta.parseall(text; filename=path)
    stmts = Tuple{Int,Any}[]
    line = 1
    for a in top.args
        if a isa LineNumberNode
            line = a.line
        else
            push!(stmts, (line, a))
        end
    end
    fns = Dict{String,Any}()
    for (k, (l0, ex)) in enumerate(stmts)
        name = _fname(ex)
        name === nothing && continue
        l1 = k < length(stmts) ? stmts[k+1][1] - 1 : nlines
        f = get!(fns, name) do
            Dict("name" => name, "line" => l0, "ranges" => UnitRange{Int}[], "doc" => "")
        end
        push!(f["ranges"], l0:l1)
        isempty(f["doc"]) && (f["doc"] = _summary(_docstring(ex)))
    end
    sort!(collect(values(fns)); by=f -> f["line"]), nlines
end
function executed(info)                              # executed lines of an LCOV tracefile: file => Set(lines)
    hits = Dict{String,Set{Int}}()
    file = ""
    for l in eachline(info)
        if startswith(l, "SF:")
            file = realpath(l[4:end])
        elseif startswith(l, "DA:")
            ln, c = parse.(Int, split(l[4:end], ',')[1:2])
            c > 0 && push!(get!(hits, file, Set{Int}()), ln)
        end
    end
    hits
end

# the code of src/ except the module file (includes, exports, the precompile workload)
files = sort([relpath(joinpath(d, f), SRC) for (d, _, fs) in walkdir(SRC) for f in fs if endswith(f, ".jl") && f != "ThetaL.jl"])
tmp = isempty(ARGS) ? mktempdir() : mkpath(ARGS[1])
if !all(X -> isfile(joinpath(tmp, X.id * ".info")), SETTINGS)
    cmd = `$(Base.julia_cmd()) -t $(Threads.nthreads()) --project=$ROOT --code-coverage=@$SRC --code-coverage=$(joinpath(tmp, "exit.info")) $(@__FILE__) $tmp`
    println("running ", length(SETTINGS), " settings under coverage …")
    run(cmd)
end
hits = Dict(X.id => executed(joinpath(tmp, X.id * ".info")) for X in SETTINGS)
fnlist = Any[]; filelist = Any[]
for f in files
    path = realpath(joinpath(SRC, f))
    fns, nlines = functions_of(path)
    push!(filelist, (path=f, dir=dirname(f) == "" ? "src" : dirname(f), lines=nlines))
    for fn in fns
        used = [X.id for X in SETTINGS if any(r -> any(in(get(hits[X.id], path, Set{Int}())), r), fn["ranges"])]
        push!(fnlist, (name=fn["name"], file=f, line=fn["line"], doc=fn["doc"], used=used))
    end
end

# ---- write the data ---------------------------------------------------------------------------------------------
js(X::AbstractString) = "\"" * escape_string(X) * "\""
js(v::AbstractVector) = "[" * join(js.(v), ",") * "]"
js(x::Real) = string(x)
js(b::Bool) = string(b)
js(nt::NamedTuple) = "{" * join((string(k) * ":" * js(v) for (k, v) in pairs(nt)), ",") * "}"
open(OUT, "w") do io
    println(io, "// Generated by docs/codemap.jl: the top-level functions of ThetaL's src/ and the settings that run them.")
    println(io, "window.CODEMAP = {")
    rev = try strip(read(`git -C $ROOT rev-parse --short HEAD`, String)) catch; "?" end
    println(io, "  generated: ", js(string(Base.Libc.strftime("%Y-%m-%d", time()), ", src at commit ", rev, " (with uncommitted changes, if any)")), ",")
    println(io, "  settings: [")
    for X in SETTINGS
        println(io, "    ", js((id=X.id, label=X.label, group=X.group, about=X.about, base=X.base, parent=X.parent, schemes=X.schemes)), ",")
    end
    println(io, "  ],\n  files: [")
    for f in filelist
        println(io, "    ", js(f), ",")
    end
    println(io, "  ],\n  functions: [")
    for f in fnlist
        println(io, "    ", js(f), ",")
    end
    println(io, "  ],\n};")
end
println("wrote ", relpath(OUT, ROOT), ": ", length(fnlist), " functions in ", length(filelist), " files, ",
        count(f -> !isempty(f.used), fnlist), " run by at least one setting")
