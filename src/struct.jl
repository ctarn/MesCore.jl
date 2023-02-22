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
