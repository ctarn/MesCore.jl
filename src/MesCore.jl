module MesCore

using Dates

fork(x::T; kwargs...) where T = T(values((; (n => getfield(x, n) for n in fieldnames(T))..., kwargs...))...)

argquery(a, lower, upper) = searchsortedfirst(a, lower):searchsortedlast(a, upper)
argquery_ε(a, x, ε) = searchsortedfirst(a, (1 - ε) * x):searchsortedlast(a, (1 + ε) * x)
argquery_δ(a, x, δ) = searchsortedfirst(a, x - δ):searchsortedlast(a, x + δ)

query(a, lower, upper) = a[argquery(a, lower, upper)]
query_ε(a, x, ε) = a[argquery_ε(a, x, ε)]
query_δ(a, x, δ) = a[argquery_δ(a, x, δ)]

error_rel(a, b) = (b - a) / max(abs(a), abs(b))
error_ppm(a, b) = error_rel(a, b) * 1e6
error_ppb(a, b) = error_rel(a, b) * 1e9
# moe: margin of error
in_moe(a, b, ε) = abs(a - b) <= max(abs(a), abs(b)) * ε

include("Ions.jl")
include("Peaks.jl")
include("MSs.jl")
include("Formulas.jl")
include("mass.jl")
include("math.jl")
include("ion.jl")
include("util.jl")
include("io/MSx.jl")
include("io/MGF.jl")
include("io/FASTA.jl")

end
