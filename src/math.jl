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
