abstract type AbstractIon end
abstract type AbstractPeak end
abstract type AbstractMS end
abstract type AbstractTandemMS <: AbstractMS end

Base.isless(a::AbstractIon, b::AbstractIon) = a.z < b.z || a.z == b.z && a.mz < b.mz

Base.isless(a::AbstractPeak, b::AbstractPeak) = a.mz < b.mz
Base.isless(peak::AbstractPeak, mz::Real) = peak.mz < mz
Base.isless(mz::Real, peak::AbstractPeak) = mz < peak.mz

Base.isless(a::AbstractMS, b::AbstractMS) = a.id < b.id

struct Ion <: AbstractIon
    m::Float64
    z::Int
    mz::Float64
    Ion(mz, z) = new(mz * z, z, mz)
    Ion(m, z, mz) = new(m, z, mz)
end

struct Peak <: AbstractPeak
    mz::Float64
    inten::Float64
end

Base.@kwdef struct MS1 <: AbstractMS
    id::Int = 0
    retention_time::Float64 = 0.0
    injection_time::Float64 = 0.0
    peaks::Vector{<:AbstractPeak} = Peak[]
end

Base.@kwdef struct MS2 <: AbstractTandemMS
    id::Int = 0
    pre::Int = 0
    retention_time::Float64 = 0.0
    injection_time::Float64 = 0.0
    activation_center::Float64 = 0.0
    isolation_width::Float64 = 0.0
    ions::Vector{<:AbstractIon} = Ion[]
    peaks::Vector{<:AbstractPeak} = Peak[]
end

dict_by_id(M::AbstractArray{<:AbstractMS}) = Dict([m.id => m for m in M])

abstract type AbstractFormula end

struct Formula <: AbstractFormula
    data::NamedTuple
end

Formula(; kwargs...) = Formula(NamedTuple(kwargs))

Base.propertynames(f::AbstractFormula, private::Bool=false) = propertynames(getfield(f, :data), private)
Base.getproperty(f::AbstractFormula, name::Symbol) = getproperty(getfield(f, :data), name)
Base.getindex(f::AbstractFormula, name::Symbol) = getproperty(f, name)

Base.parse(::Type{Formula}, s::AbstractString) = begin
    es = map(e -> split(e, '('), split(s, ')'; keepempty=false))
    es = map(((e, n),) -> (e=Symbol(e), n=parse(Float64, n)), es)
    d = Dict{Symbol, Float64}()
    foreach(x -> d[x.e] = get(d, x.e, 0) + x.n, es)
    return Formula(; d...)
end

calc_mass(f::Formula, tab::NamedTuple) = sum(f[n] * tab[n] for n in propertynames(f); init=0.0)
