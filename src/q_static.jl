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
