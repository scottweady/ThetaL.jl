# Terminal plots of θ–L curves with DotPlots (loaded with `using ThetaL, DotPlots`): a curve or a vector of curves is
# one frame, an `Evolution` one frame per saved time, colored by time. The keywords are DotPlots' (`animate`,
# `every`, `width`, …).
module ThetaLDotPlotsExt

using ThetaL: ThetaLState, Evolution, curve_points
import DotPlots: dotplot

const _Curves = Union{ThetaLState,AbstractVector{ThetaLState}}

_frame(X::ThetaLState) = [curve_points(X)]
_frame(Xs::AbstractVector{ThetaLState}) = [curve_points(X) for X in Xs]

dotplot(X::_Curves; kw...) = dotplot(stdout, X; kw...)
dotplot(io::IO, X::_Curves; kw...) = dotplot(io, [_frame(X)]; kw...)
dotplot(sol::Evolution; kw...) = dotplot(stdout, sol; kw...)
dotplot(io::IO, sol::Evolution; kw...) = dotplot(io, [_frame(X) for X in sol.curves]; t=sol.t, kw...)

end
