abstract type AbstractPeak end

Base.isless(a::AbstractPeak, b::AbstractPeak) = a.mz < b.mz
Base.isless(peak::AbstractPeak, mz::Real) = peak.mz < mz
Base.isless(mz::Real, peak::AbstractPeak) = mz < peak.mz

struct Peak <: AbstractPeak
    mz::Float64
    inten::Float64
end
