# These core functions are designed to work with Vector{T}, StaticVector{T}, and with
# Tuple{T, T, T, T}. They return scalars (i.e., tuple)

# Core quaternion multiplication function, implementing p ⊗ q
function q_multiplyCore(p, q)
    ps, px, py, pz = p
    qs, qx, qy, qz = q
    os = ps * qs - px * qx - py * qy - pz * qz
    ox = px * qs + ps * qx - pz * qy + py * qz
    oy = py * qs + pz * qx + ps * qy - px * qz
    oz = pz * qs - py * qx + px * qy + ps * qz
    return os, ox, oy, oz
end

# R_AB from q_AB
function q_toDcmCore(q)
    qs, qx, qy, qz = q
    x2, y2, z2 = qx + qx, qy + qy, qz + qz
    sx, sy, sz = qs * x2, qs * y2, qs * z2
    xx, xy, xz = qx * x2, qx * y2, qx * z2
    yy, yz, zz = qy * y2, qy * z2, qz * z2

    R11 = 1 - (yy + zz)
    R12 = xy - sz
    R13 = xz + sy

    R21 = xy + sz
    R22 = 1 - (xx + zz)
    R23 = yz - sx

    R31 = xz - sy
    R32 = yz + sx
    R33 = 1 - (xx + yy)
    return R11, R21, R31, R12, R22, R32, R13, R23, R33
end

# R_AB to q_AB
function q_fromDcmCore(R)
    r11, r21, r31, r12, r22, r32, r13, r23, r33 = R
    vmax = 1 + r11 - r22 - r33
    v2 = 1 - r11 + r22 - r33
    v3 = 1 - r11 - r22 + r33
    v4 = 1 + r11 + r22 + r33

    idx = 1
    if v2 > vmax
        idx = 2
        vmax = v2
    end
    if v3 > vmax
        idx = 3
        vmax = v3
    end
    if v4 > vmax
        idx = 4
        vmax = v4
    end

    qq = sqrt(abs(vmax)) / 2
    f = 0.25 / qq

    idx == 1 && return f * (r32 - r23), qq, f * (r12 + r21), f * (r31 + r13)
    idx == 2 && return f * (r13 - r31), f * (r12 + r21), qq, f * (r23 + r32)
    idx == 3 && return f * (r21 - r12), f * (r31 + r13), f * (r23 + r32), qq
    return -qq, f * (r23 - r32), f * (r31 - r13), f * (r12 - r21)
end

# v_A from q_AB, v_B
function q_transformVectorCore(q_AB, v_B)
    qs, qx, qy, qz = q_AB
    x, y, z = v_B
    cx = qy * z - qz * y
    cy = qz * x - qx * z
    cz = qx * y - qy * x

    c2x = qy * cz - qz * cy
    c2y = qz * cx - qx * cz
    c2z = qx * cy - qy * cx

    xo = x + 2 * (c2x + qs * cx)
    yo = y + 2 * (c2y + qs * cy)
    zo = z + 2 * (c2z + qs * cz)
    return xo, yo, zo
end

@inline function q_transposeCore(q)
    qs, qx, qy, qz = q
    return qs, -qx, -qy, -qz
end

function q_inverseCore(q)
    qs, qx, qy, qz = q
    nq2 = qs * qs + qx * qx + qy * qy + qz * qz
    return qs / nq2, -qx / nq2, -qy / nq2, -qz / nq2
end

# dq_AB from q_AB, ωAB_B
function q_derivativeCore(q, w)
    qs, qx, qy, qz = q
    wx, wy, wz = w
    dqs = (-qx * wx - qy * wy - qz * wz) / 2
    dqx = (qs * wx - qz * wy + qy * wz) / 2
    dqy = (qz * wx + qs * wy - qx * wz) / 2
    dqz = (-qy * wx + qx * wy + qs * wz) / 2
    return dqs, dqx, dqy, dqz
end

@inline function q_toRvCore(q)
    qs, qx, qy, qz = q
    nqv = sqrt(qx*qx + qy*qy + qz*qz)
    iszero(nqv) && return zero(qx), zero(qy), zero(qz)
    f = 2*atan(nqv, qs) / nqv
    return f*qx, f*qy, f*qz
end
