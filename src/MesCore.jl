module MesCore

fork(x::T; kwargs...) where T = T(values((; (n => getfield(x, n) for n in fieldnames(T))..., kwargs...))...)

mapkey(f, D::Dict) = zip(f.(keys(D)), values(D)) |> Dict
mapvalue(f, D::Dict) = zip(keys(D), f.(values(D))) |> Dict
mapkey(f) = Base.Fix1(mapkey, f)
mapvalue(f) = Base.Fix1(mapvalue, f)

argquery(a, lower, upper) = searchsortedfirst(a, lower):searchsortedlast(a, upper)
argquery_δ(a, x, δ) = argquery(a, x - δ, x + δ)
argquery_ε(a, x, ε) = argquery_δ(a, x, ε * x)

query(a, lower, upper) = a[argquery(a, lower, upper)]
query_δ(a, x, δ) = query(a, x - δ, x + δ)
query_ε(a, x, ε) = query_δ(a, x, ε * x)

argquery_near(a, x) = begin
    i = searchsortedfirst(a, x)
    if i > lastindex(a)
        return lastindex(a)
    elseif i == firstindex(a)
        return firstindex(a)
    else
        return i - (abs(a[i] - x) >= abs(a[i-1] - x))
    end
end

query_near(a, x) = a[argquery_near(a, x)]

include("Consts.jl")
include("Ions.jl")
include("Peaks.jl")
include("MSs.jl")
include("Formulas.jl")
include("math.jl")
include("util.jl")
include("io/MSx.jl")
include("io/MGF.jl")
include("io/FASTA.jl")

end
