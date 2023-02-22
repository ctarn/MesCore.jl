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

calc_centroid(xs, ws) = begin
    if length(xs) > 1
        ss = (ws[begin:end-1] .+ ws[begin+1:end]) .* (xs[begin+1:end] .- xs[begin:end-1])
        s = sum(ss) / 2
        x = sum(ss .* (xs[begin:end-1] .+ xs[begin+1:end])) / 4
        return x / s, s / (xs[end] - xs[begin])
    elseif length(xs) == 1
        return xs[begin], ws[begin]
    else
        return zero(eltype(xs)), zero(eltype(ws))
    end
end

log_softer(s=1) = x -> sign(x) * s * (log(abs(x) + s) - log(s))

include("Ions.jl")
include("Peaks.jl")
include("MSs.jl")
include("Formulas.jl")
include("mass.jl")
include("io/MSx.jl")
include("io/MGF.jl")
include("io/FASTA.jl")
include("ion.jl")

match_path(path, ext="") = begin
    paths = readdir(dirname(path); join=true)
    return filter(p -> startswith(basename(p), basename(path)) && endswith(p, ext), paths)
end

parse_range(::Type{T}, s::AbstractString) where T = begin
    vs = map(x -> parse(T,  x),  split(s, ':'))
    if length(vs) <= 2 return range(vs[begin], vs[end])
    elseif length(vs) == 3 return range(vs[begin], vs[end]; step=vs[begin+1])
    else return nothing
    end
end

date_mark(date=today()) = replace(string(date), "-"=>"")

open_url(url) = begin
    if Sys.isapple()
        run(`open $(url)`)
    elseif Sys.iswindows()
        run(`cmd /k start $(url) "&&" exit`)
    else
        run(`open $(url)`)
    end
end

end
