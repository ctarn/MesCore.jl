const mₚ = 1.007276466621

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
