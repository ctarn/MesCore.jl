using Statistics

error_rel(a, b) = (b - a) / max(abs(a), abs(b))
error_ppm(a, b) = error_rel(a, b) * 1e6
error_ppb(a, b) = error_rel(a, b) * 1e9

# moe: margin of error
in_moe(a, b, ε) = abs(a - b) <= max(abs(a), abs(b)) * ε

m_to_mz(m, z) = (m + mₚ * z) / abs(z)
mz_to_m(mz, z) = mz * abs(z) - mₚ * z

mh_to_mz(mh, z) = begin
    if z > 0 return (mh + mₚ * (z - 1)) / z end
    if z < 0 return (mh - mₚ * (-z - 1)) / -z end
    return 0.0
end

mz_to_mh(mz, z) = begin
    if z > 0 return (mz - mₚ) * z + mₚ end
    if z < 0 return (mz + mₚ) * -z - mₚ end
    return 0.0
end

calc_centroid(xs, ws) = begin
    @assert length(xs) == length(ws) "`xs` and `ws` are not of the same length"
    @assert !isempty(xs) "`xs` and `ws` are both empty"
    if xs[begin] == xs[end]
        return xs[begin], mean(ws)
    else
        ss2 = (ws[begin:end-1] .+ ws[begin+1:end]) .* (xs[begin+1:end] .- xs[begin:end-1])
        s2 = sum(ss2)
        x2 = sum(ss2 .* (xs[begin:end-1] .+ xs[begin+1:end])) / 2
        return x2 / s2, s2 / 2 / (xs[end] - xs[begin])
    end
end

log_softer(s=1) = x -> copysign(s * (log(abs(x) + s) - log(s)), x)
exp_softer(s=1) = x -> x * exp(-abs(x)/s)
