error_rel(a, b) = (b - a) / max(abs(a), abs(b))
error_ppm(a, b) = error_rel(a, b) * 1e6
error_ppb(a, b) = error_rel(a, b) * 1e9

# moe: margin of error
in_moe(a, b, ε) = abs(a - b) <= max(abs(a), abs(b)) * ε

m_to_mz(m, z) = (m + mₚ * z) / abs(z)
mz_to_m(mz, z) = mz * abs(z) - mₚ * z

mh_to_mz(mh, z) = begin
    if z > 0
        (mh + mₚ * (z - 1)) / z
    elseif z < 0
        (mh - mₚ * (-z - 1)) / -z
    else
        0.0
    end
end

mz_to_mh(mz, z) = begin
    if z > 0
        mz * z - mₚ * (z - 1)
    elseif z < 0
        mz * -z + mₚ * (-z - 1)
    else
        0.0
    end
end

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
