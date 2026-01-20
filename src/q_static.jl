@inline function q_toDcm(q::SVector{4,T}) where {T}
    s, x, y, z = q
    x2, y2, z2 = x + x, y + y, z + z
    sx, sy, sz = s*x2, s*y2, s*z2
    xx, xy, xz = x*x2, x*y2, x*z2
    yy, yz, zz = y*y2, y*z2, z*z2

    return @SMatrix [
        (1 - (yy + zz)) (xy - sz) (xz + sy);
        (xy + sz) (1 - (xx + zz)) (yz - sx);
        (xz - sy) (yz + sx) (1 - (xx + yy))
    ]
end

@inline function q_fromDcm(R_BA::SMatrix{3,3,T}) where {T}
    r11, r21, r31 = R_BA[1, 1], R_BA[1, 2], R_BA[1, 3]
    r12, r22, r32 = R_BA[2, 1], R_BA[2, 2], R_BA[2, 3]
    r13, r23, r33 = R_BA[3, 1], R_BA[3, 2], R_BA[3, 3]

    vmax = 1 + r11 - r22 - r33
    v2 = 1 - r11 + r22 - r33
    v3 = 1 - r11 - r22 + r33
    v4 = 1 + r11 + r22 + r33

    idx = 1
    if v2 > vmax
        ;
        idx = 2;
        vmax = v2;
    end
    if v3 > vmax
        ;
        idx = 3;
        vmax = v3;
    end
    if v4 > vmax
        ;
        idx = 4;
        vmax = v4;
    end

    qx = 0.5*sqrt(abs(vmax))
    f = 0.25/qx

    if idx == 1
        return @SVector [f*(r23 - r32); qx; f*(r12 + r21); f*(r31 + r13)]
    elseif idx == 2
        return @SVector [f*(r31 - r13); f*(r12 + r21); qx; f*(r23 + r32)]
    elseif idx == 3
        return @SVector [f*(r12 - r21); f*(r31 + r13); f*(r23 + r32); qx]
    end
    return @SVector [qx; f*(r23 - r32); f*(r31 - r13); f*(r12 - r21)]
end

@inline function q_transpose(q::SVector{4,T}) where {T}
    return @SVector [q[1]; -q[2]; -q[3]; -q[4]]
end

@inline function q_derivative(q_AB::SVector{4,T}, ωAB_B::SVector{3,T}) where {T}
    ps, px, py, pz = q_AB
    qx, qy, qz = ωAB_B

    return @SVector [
        (- px*qx - py*qy - pz*qz)/2;
        (+ ps*qx - pz*qy + py*qz)/2;
        (+ pz*qx + ps*qy - px*qz)/2;
        (- py*qx + px*qy + ps*qz)/2
    ]
end

@inline function q_tox(q_AB::SVector{4,T}) where {T}
    qs, qx, qy, qz = q_AB

    c2x = -qy*qy - qz*qz
    c2y = qx*qy
    c2z = qx*qz

    return @SVector [
        1.0 + 2*c2x;
        2(c2y + qs*qz);
        2(c2z - qs*qy)
    ]
end

@inline function q_toy!(q_AB::SVector{4,T}) where {T}
    qs, qx, qy, qz = q_AB

    c2x = qy*qx
    c2y = -qz*qz - qx*qx
    c2z = qy*qz

    return @SVector [
        2(c2x - qs*qz);
        1.0 + 2(c2y);
        2(c2z + qs*qx)
    ]
end

@inline function q_toz!(q_AB::SVector{4,T}) where {T}
    qs, qx, qy, qz = q_AB

    c2x = qz*qx
    c2y = qz*qy
    c2z = -qx*qx - qy*qy

    return @SVector [
        2(c2x + qs*qy);
        2(c2y - qs*qx);
        1.0 + 2(c2z)
    ]
end

@inline function q_attitudeError(qNominal::SVector{4,T}, q::SVector{4,T}) where {T}
    ps, px, py, pz = q
    qs, qx, qy, qz = qNominal

    imax = findmax(abs, qNominal)[2]
    sgn = sign(qNominal[imax]) == sign(q[imax]) ? 2.0 : -2.0

    return @SVector [
        (-px*qs + ps*qx + pz*qy - py*qz)*sgn;
        (-py*qs - pz*qx + ps*qy + px*qz)*sgn;
        (-pz*qs + py*qx - px*qy + ps*qz)*sgn
    ]
end
