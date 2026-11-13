# ThetaL regression tests (test_*.jl here), no solve: states of curves (test_curves.jl) and θ–L steps under
# prescribed velocities and symbols (test_stepping.jl). Each file runs in its own module; shared test curves are in
# common/.
# From this package's directory:
#   julia -t 8 --project -e 'using Pkg; Pkg.test()'                                   (all files)
#   julia -t 8 --project -e 'using Pkg; Pkg.test(test_args=["stepping"])'             (only test_stepping.jl)
using Test, Printf

files = sort(filter(f -> startswith(f, "test_") && endswith(f, ".jl"), readdir(@__DIR__)))
isempty(ARGS) || (files = filter(f -> f[6:end-3] in ARGS, files))

@testset "ThetaL" begin
    for f in files
        t = @elapsed @testset verbose=true "$f" begin
            m = Core.eval(Main, :(module $(Symbol(f[1:end-3])) end))
            Base.include(m, joinpath(@__DIR__, f))
        end
        @printf("  %-30s %7.1f s\n", f, t)
    end
end
