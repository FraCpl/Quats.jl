"""
    q_AC = q_multiply(q_AB, q_BC)

Mutliply the two input quaternions as follows:
```math
q_{AC} = q_{AB} ⊗ q_{BC}
```
"""
@inline function q_multiply(q_AB, q_BC)
    q_AC = similar(q_AB)
    q_multiply!(q_AC, q_AB, q_BC)
    return q_AC
end

# Baseline in-place multiplication function
@inline function q_multiply!(q_AC, q_AB, q_BC)
    ps, px, py, pz = q_AB
    qs, qx, qy, qz = q_BC
    q_multiplyCore!(q_AC, ps, px, py, pz, qs, qx, qy, qz)
    return nothing
end

# Auiliary multiplication functions: these are included to allow high-speed in-place
# multiplication operations without the need of allocating/computing a quaternion transpose

# First term of multiplication is the transpose of the baseline one
@inline function q_multiplyT1!(q_AC, q_BA, q_BC)
    ps, px, py, pz = q_BA
    qs, qx, qy, qz = q_BC
    q_multiplyCore!(q_AC, ps, -px, -py, -pz, qs, qx, qy, qz)
    return nothing
end

# Second term of multiplication is the transpose of the baseline one
@inline function q_multiplyT2!(q_AC, q_AB, q_CB)
    ps, px, py, pz = q_AB
    qs, qx, qy, qz = q_CB
    q_multiplyCore!(q_AC, ps, px, py, pz, qs, -qx, -qy, -qz)
    return nothing
end

# Both first and second terms of multiplication are the transpose of the baseline ones
@inline function q_multiplyT12!(q_AC, q_BA, q_CB)
    ps, px, py, pz = q_BA
    qs, qx, qy, qz = q_CB
    q_multiplyCore!(q_AC, ps, -px, -py, -pz, qs, -qx, -qy, -qz)
    return nothing
end

# q1 <- q1 * q2
@inline function q_multiply!(q1, q2)
    ps, px, py, pz = q1
    qs, qx, qy, qz = q2
    q_multiplyCore!(q1, ps, px, py, pz, qs, qx, qy, qz)
    return nothing
end

# Core quaternion multiplication function
@inline function q_multiplyCore!(qOut, ps, px, py, pz, qs, qx, qy, qz)
    # p ⊗ q
    qOut[1] = ps*qs - px*qx - py*qy - pz*qz
    qOut[2] = px*qs + ps*qx - pz*qy + py*qz
    qOut[3] = py*qs + pz*qx + ps*qy - px*qz
    qOut[4] = pz*qs - py*qx + px*qy + ps*qz
    return nothing
end

"""
    q = q_multiplyn(q1, q2, q3, ...)

```math
    q = q₁ ⊗ q₂ ⊗ q₃ ⊗ ...
```
"""
@inline function q_multiplyn(q...)
    qOut = copy(q[1])

    for i in 2:lastindex(q)
        qOut = q_multiply(qOut, q[i])
    end

    return qOut
end

"""
    q_build(qs,qv)

Build a quaternion from its scalar and vectorial components.
"""
@inline q_build(qs, qv) = [qs; qv]

"""
    R_AB = q_toDcm(q_AB)

Translate the input unitary quaternion into a transformation matrix.
```math
R_{AB}(q_{AB}) = I + 2qₛ[qᵥ×] + 2[qᵥ×]²
```
"""
@inline function q_toDcm(q)
    R_BA = Matrix{Float64}(undef, 3, 3)
    q_toDcm!(R_BA, q)
    return R_BA
end

@views @inline function q_toDcm!(R, q)     # R_BA from q_BA
    # qx = crossMat(q[2:4])
    # return I + 2.0*(qx*qx + q[1].*qx)
    s, x, y, z = q
    x2, y2, z2 = x + x, y + y, z + z
    sx, sy, sz = s*x2, s*y2, s*z2
    xx, xy, xz = x*x2, x*y2, x*z2
    yy, yz, zz = y*y2, y*z2, z*z2

    R[1, 1] = 1.0 - (yy + zz)
    R[1, 2] = xy - sz
    R[1, 3] = xz + sy

    R[2, 1] = xy + sz
    R[2, 2] = 1.0 - (xx + zz)
    R[2, 3] = yz - sx

    R[3, 1] = xz - sy
    R[3, 2] = yz + sx
    R[3, 3] = 1.0 - (xx + yy)

    return nothing
end

"""
    q_AB = q_fromDcm(R_AB)

Translate the input rotation matrix into a unitary quaternion.
"""
@inline function q_fromDcm(R_AB)
    q_AB = Vector{eltype(R_AB)}(undef, 4)
    q_fromDcm!(q_AB, R_AB)
    return q_AB
end

@inline function q_fromDcm!(q_AB, R_AB)
    r11, r21, r31 = R_AB[1, 1], R_AB[1, 2], R_AB[1, 3]
    r12, r22, r32 = R_AB[2, 1], R_AB[2, 2], R_AB[2, 3]
    r13, r23, r33 = R_AB[3, 1], R_AB[3, 2], R_AB[3, 3]
    q_fromDcmCore!(q_AB, r11, r12, r13, r21, r22, r23, r31, r32, r33)
    return nothing
end

@inline function q_fromDcmCore!(q, r11, r12, r13, r21, r22, r23, r31, r32, r33)
    # dcm11 = R_BA[1, 1]; dcm12 = R_BA[2, 1]; dcm13 = R_BA[3, 1];
    # dcm21 = R_BA[1, 2]; dcm22 = R_BA[2, 2]; dcm23 = R_BA[3, 2];
    # dcm31 = R_BA[1, 3]; dcm32 = R_BA[2, 3]; dcm33 = R_BA[3, 3];

    # vv = 1.0 .+ [+dcm11-dcm22-dcm33; -dcm11+dcm22-dcm33; -dcm11-dcm22+dcm33; +dcm11+dcm22+dcm33]
    # idx = argmax(vv);
    # qx = 0.5*sqrt(abs(vv[idx]))
    # f  = 0.25/qx;

    # if idx == 1
    #     return [f*(dcm23 - dcm32); qx; f*(dcm12 + dcm21);  f*(dcm31 + dcm13)]
    # elseif idx == 2
    #     return [f*(dcm31 - dcm13); f*(dcm12 + dcm21); qx; f*(dcm23 + dcm32)]
    # elseif idx == 3
    #     return [f*(dcm12 - dcm21); f*(dcm31 + dcm13); f*(dcm23 + dcm32); qx]
    # end
    # return [qx; f*(dcm23 - dcm32); f*(dcm31 - dcm13); f*(dcm12 - dcm21)]

    # r11, r21, r31 = R_BA[1, 1], R_BA[1, 2], R_BA[1, 3]
    # r12, r22, r32 = R_BA[2, 1], R_BA[2, 2], R_BA[2, 3]
    # r13, r23, r33 = R_BA[3, 1], R_BA[3, 2], R_BA[3, 3]

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

    qx = 0.5*sqrt(abs(vmax))
    f = 0.25/qx

    if idx == 1
        q[1] = f*(r23 - r32)
        q[2] = qx
        q[3] = f*(r12 + r21)
        q[4] = f*(r31 + r13)
    elseif idx == 2
        q[1] = f*(r31 - r13)
        q[2] = f*(r12 + r21)
        q[3] = qx
        q[4] = f*(r23 + r32)
    elseif idx == 3
        q[1] = f*(r12 - r21)
        q[2] = f*(r31 + r13)
        q[3] = f*(r23 + r32)
        q[4] = qx
    else
        q[1] = qx
        q[2] = f*(r23 - r32)
        q[3] = f*(r31 - r13)
        q[4] = f*(r12 - r21)
    end
    return nothing
end

"""
    q_AB = q_fromAxes(xB_A, yB_A, zB_A)

Compute the attitude quaternion given as input the axes of a reference frame.
"""
# @inline q_fromAxes(xB_A, yB_A, zB_A) = q_fromDcm(dcm_fromAxes(xB_A, yB_A, zB_A))

@inline function q_fromAxes(xB_A, yB_A, zB_A)
    q_AB = Vector{eltype(yB_A)}(undef, 4)
    q_fromAxes!(q_AB, xB_A, yB_A, zB_A)
    return q_AB
end

@inline function q_fromAxes!(q_AB, xB_A, yB_A, zB_A)
    r11, r12, r13 = xB_A
    r21, r22, r23 = yB_A
    r31, r32, r33 = zB_A
    q_fromDcmCore!(q_AB, r11, r12, r13, r21, r22, r23, r31, r32, r33)
    return nothing
end

"""
    q_random()

Generate a random unitary quaternion.
"""
@inline q_random() = normalize(randn(4))

"""
    q_identity()

Get identity quaternion, with scalar component equal to 1 and
vector components equal to zero.
"""
@inline q_identity() = [1.0; 0.0; 0.0; 0.0]

"""
    v_A = q_transformVector(q_AB, v_B)

Project the vector v from frame B into frame A.
"""
@inline function q_transformVector(q_AB, v_B)
    v_A = similar(v_B)
    q_transformVector!(v_A, q_AB, v_B)
    return v_A
end

"""
    q_transformVector!(v_A, q_AB, v_B)

Project the vector v from frame B into frame A.
"""
@inline function q_transformVector!(v_A, q_AB, v_B)
    # qxv = q_AB[2:4] × v_B
    # return v_B + 2.0*(q_AB[2:4] × qxv + q_AB[1].*qxv) # v_A
    qs, qx, qy, qz = q_AB
    q_transformVectorCore!(v_A, qs, qx, qy, qz, v_B)
    return nothing
end

"""
    q_transformVectorT!(v_A, q_BA, v_B)

Project the vector v from frame B into frame A, using q_BA.
"""
@inline function q_transformVectorT!(v_A, q_BA, v_B)
    qs, qx, qy, qz = q_BA
    q_transformVectorCore!(v_A, qs, -qx, -qy, -qz, v_B)
    return nothing
end

# q = q_AB
@inline function q_transformVectorCore!(Vout, qs, qx, qy, qz, Vin)
    x, y, z = Vin

    cx = qy*z - qz*y
    cy = qz*x - qx*z
    cz = qx*y - qy*x

    c2x = qy*cz - qz*cy
    c2y = qz*cx - qx*cz
    c2z = qx*cy - qy*cx

    Vout[1] = x + 2(c2x + qs*cx)
    Vout[2] = y + 2(c2y + qs*cy)
    Vout[3] = z + 2(c2z + qs*cz)
    return nothing
end

"""
    q' = q_transpose(q)

Transpose the input quaternion.
"""
@inline function q_transpose!(q)
    @inbounds for i in 2:4
        ;
        q[i] = -q[i];
    end
    return nothing
end
@inline q_transpose(q) = [q[1]; -q[2]; -q[3]; -q[4]]
@inline function q_transpose!(qt, q)
    qt[1] = q[1]
    @inbounds for i in 2:4
        ;
        qt[i] = -q[i];
    end
    return nothing
end

"""
    q̇_AB = q_derivative(q_AB, ωAB_B)

Compute the time derivative of a unitary quaternion, given the corresponding
angular velocity vector.

Mathematically, this function performs the following operation:
```math
q̇_{AB} = \\frac{1}{2} q_{AB} ⊗ [0; ω^B_{AB}]
```
where ```ωAB_B``` represents the angular velocity of frame ``B`` with respect to
frame ``A``, projected into frame ``B``.
"""
@inline function q_derivative(q_AB, ωAB_B)
    dq_AB = similar(q_AB)
    q_derivative!(dq_AB, q_AB, ωAB_B)
    return dq_AB
end

"""
    q_derivative!(q̇_AB, q_AB, ωAB_B)

Compute the time derivative of a unitary quaternion, given the corresponding
angular velocity vector.

Mathematically, this function performs the following operation:
```math
q̇_{AB} = \\frac{1}{2} q_{AB} ⊗ [0; ω^B_{AB}]
```
where ```ωAB_B``` represents the angular velocity of frame ``B`` with respect to
frame ``A``, projected into frame ``B``.
"""
@inline function q_derivative!(dq_AB, q_AB, ωAB_B)
    # q_multiply(q_AB, [0.0; 0.5.*ωAB_B])  # dq_BA

    ps, px, py, pz = q_AB
    qx, qy, qz = ωAB_B

    dq_AB[1] = (-px * qx - py * qy - pz * qz) / 2
    dq_AB[2] = (ps * qx - pz * qy + py * qz) / 2
    dq_AB[3] = (pz * qx + ps * qy - px * qz) / 2
    dq_AB[4] = (-py * qx + px * qy + ps * qz) / 2

    return nothing
end

"""
    q_AB = q_fromAxisAngle(u, θ_AB)

Compute the unitary quaternion given as input an axis-angle representation.
"""
q_fromAxisAngle(u, θ) = iszero(u) ? q_identity() : [cos(0.5θ); sin(0.5θ)*normalize(u)]

@inline function q_fromAxisAngle(idx::Int, θ)
    sθ, cθ = sincos(θ/2)
    q = [cθ; 0.0; 0.0; 0.0]
    q[idx + 1] = sθ
    return q
end

function q_fromAxisAngle!(q, idx::Int, θ)
    sθ, cθ = sincos(θ / 2)
    q .= 0
    q[1] = cθ
    q[idx + 1] = sθ
    return nothing
end

function q_fromAxisAngle!(q, u, θ)
    uNormSq = 0.0
    @inbounds for i in 1:3
        uNormSq += u[i] * u[i]
    end

    if uNormSq == 0
        q .= 0.0
        q[1] = 1.0
        return nothing
    end

    st, q[1] = sincos(θ / 2)
    st = st / sqrt(uNormSq)
    @inbounds for i in 1:3
        q[i + 1] = st * u[i]
    end
    return nothing
end

@inline function q_toAxisAngle(q)
    qs, qx, qy, qz = q
    nqv = sqrt(qx*qx + qy*qy + qz*qz)
    return [qx/nqv; qy/nqv; qz/nqv], 2atan(nqv, qs)
end

@inline function q_toAxes(q_AB)
    R = q_toDcm(q_AB)
    return R[1, :], R[2, :], R[3, :]    # xB_A, yB_A, zB_A
end


"""
    q⁻¹ = q_inverse(q)

Compute the inverse of the input quaternion.
"""
@inline function q_inverse(q)
    qs, qx, qy, qz = q
    nq2 = qs*qs + qx*qx + qy*qy + qz*qz
    return [qs/nq2; -qx/nq2; -qy/nq2; -qz/nq2]
end

@inline function q_toRv(q)
    qs, qx, qy, qz = q
    nqv = sqrt(qx*qx + qy*qy + qz*qz)
    if nqv < 1e-10
        ;
        return zeros(3);
    end
    nqv = 2atan(nqv, qs)/nqv
    return [qx*nqv; qy*nqv; qz*nqv]
end

@inline q_fromRv(ϕ) = q_fromAxisAngle(ϕ, norm(ϕ))

@inline q_fromRv!(q, ϕ) = q_fromAxisAngle!(q, ϕ, norm(ϕ))

"""
This uses rotation vector to compute the average angular rate between two
sampled quaternion values, knowing that:
q[k] = q[k-1] o dq
and
angRate(t) = dtheta/dt for t in [t[k-1],t[k]]
This means that angRate is the equivalent constant average rate that
rotates q[k-1] into q[k].
"""
@inline @views function q_rate(t, q_AB)
    dt = diff(t)
    dq_AB = q_multiply.(q_transpose.(q_AB[1:(end - 1)]), q_AB[2:end])
    return [q_toRv.(dq_AB) ./ dt; [zeros(3)]]    # angRateAB_B
end

# Δt = t2 - t1
# q1_AB = q_AB[t1]
# q2_AB = q_AB[t2]
@inline @views function q_rate(Δt, q1_AB, q2_AB)
    dq_AB = q_multiply(q_transpose(q1_AB), q2_AB)
    return q_toRv(dq_AB)/Δt
end

# https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
# 3(yaw)-2(pitch)-1(roll) sequence
# CAUTION: Only [3, 2, 1] and [1, 2, 3] sequencs are implemented for now
function q_toEuler(q, sequence=[3, 2, 1])
    qs, qx, qy, qz = q

    if sequence == [3, 2, 1]
        # roll (x-axis rotation)
        t1 = 2(qs*qx + qy*qz)
        t2 = 1 - 2(qx*qx + qy*qy)
        roll = atan(t1, t2)

        # pitch (y-axis rotation)
        t3 = 2(qs*qy - qz*qx)
        t3 = max(-1.0, min(t3, 1.0))
        pitch = asin(t3)

        # yaw (z-axis rotation)
        t4 = 2(qs*qz + qx*qy)
        t5 = 1 - 2(qy*qy + qz*qz)
        yaw = atan(t4, t5)

        return [yaw; pitch; roll]

    elseif sequence == [1, 2, 3]
        sinPitch = 2(qx*qz + qy*qs)
        mask = sinPitch ≥ 1 - 10*eps() || sinPitch ≤ -1 + 10*eps()
        sinPitch = max(-1.0, min(sinPitch, 1.0))
        pitch = asin(sinPitch)
        if mask
            roll = 0
            yaw = sign(sinPitch)*2*atan(qx, qs)
        else
            roll = atan(-2(qy*qz - qx*qs), qs^2 - qx^2 - qy^2 + qz^2);
            yaw = atan(-2(qx*qy - qz*qs), qs^2 + qx^2 - qy^2 - qz^2)
        end
        return [roll; pitch; yaw]
    end
end

@inline function q_fromEuler(θ, sequence=[3, 2, 1])
    q = q_fromAxisAngle(sequence[1], θ[1])
    for i in 2:lastindex(θ)
        q .= q_multiply(q, q_fromAxisAngle(sequence[i], θ[i]))
    end
    return q
end

"""
    qe = q_exp(q)

Compute the exponential of the input quaternion.
"""
@inline function q_exp(q)
    qs, qx, qy, qz = q
    qvn = sqrt(qx*qx + qy*qy + qz*qz)
    qvnn = qvn == 0 ? 1.0 : qvn
    k = exp(qs)
    qvnn = k*sin(qvn)/qvnn
    return [k*cos(qvn); qx*qvnn; qy*qvnn; qy*qvnn]
end

"""
    ql = q_log(q)

Compute the logarithm of the input quaternion.
"""
@inline function q_log(q)
    qs, qx, qy, qz = q
    vNorm = sqrt(qx*qx + qy*qy + qz*qz)
    if vNorm == 0.0
        vNorm = 1.0
    end
    nq = norm(q)
    vNorm = acos(qs/nq)/vNorm
    return [log(nq); qx*vNorm; qy*vNorm; qz*vNorm]
end

@inline q_power(q, n) = q_exp(n .* q_log(q))

"""
    q = q_slerp(q0, q1, τ)

Compute the spherical linear interpolation between two quaternions at τ, where q0 = q[τ=0], q1 = q[τ=1], and q = q[τ].
"""
@inline q_slerp(q0, q1, τ) = q_multiply(q0, q_power(q_multiply(q_transpose(q0), q1), τ))

@inline function q_interp(t, q, ti)
    id0 = findlast(ti .≥ t)
    id1 = min(id0+1, length(t))
    return q_slerp(q[id0], q[id1], (ti - t[id0])/(t[id1] - t[id0]))
end

# e.g., qNominal = qDes, q = qEst
@inline function q_attitudeError(qNominal, q)
    err = Vector{eltype(q)}(undef, 3)
    q_attitudeError!(err, qNominal, q)
    return err
end

@inline function q_attitudeError!(err, qNominal, q)
    ps, px, py, pz = q
    qs, qx, qy, qz = qNominal

    imax = findmax(abs, qNominal)[2]
    sgn = sign(qNominal[imax]) == sign(q[imax]) ? 2.0 : -2.0

    err[1] = (-px*qs + ps*qx + pz*qy - py*qz)*sgn
    err[2] = (-py*qs - pz*qx + ps*qy + px*qz)*sgn
    err[3] = (-pz*qs + py*qx - px*qy + ps*qz)*sgn

    return nothing  # 2q_multiply(q_transpose(sgn*q), qNominal)[2:4]
end

"""
    xB_A = q_tox(q_AB)

Compute the x axis of frame B projected in frame A.
"""
@inline function q_tox(q_AB)
    out = Vector{eltype(q_AB)}(undef, 3)
    q_tox!(out, q_AB)
    return out
end

"""
    yB_A = q_toy(q_AB)

Compute the y axis of frame B projected in frame A.
"""
@inline function q_toy(q_AB)
    out = Vector{eltype(q_AB)}(undef, 3)
    q_toy!(out, q_AB)
    return out
end

"""
    zB_A = q_toz(q_AB)

Compute the z axis of frame B projected in frame A.
"""
@inline function q_toz(q_AB)
    out = Vector{eltype(q_AB)}(undef, 3)
    q_toz!(out, q_AB)
    return out
end

"""
    q_tox!(xB_A, q_AB)

Compute the x axis of frame B projected in frame A.
"""
@inline function q_tox!(xB_A, q_AB)
    qs, qx, qy, qz = q_AB

    c2x = -qy*qy - qz*qz
    c2y = qx*qy
    c2z = qx*qz

    xB_A[1] = 1.0 + 2*c2x
    xB_A[2] = 2(c2y + qs*qz)
    xB_A[3] = 2(c2z - qs*qy)

    return nothing
end

"""
    q_toy!(yB_A, q_AB)

Compute the y axis of frame B projected in frame A.
"""
@inline function q_toy!(yB_A, q_AB)
    qs, qx, qy, qz = q_AB

    c2x = qy*qx
    c2y = -qz*qz - qx*qx
    c2z = qy*qz

    yB_A[1] = 2(c2x - qs*qz)
    yB_A[2] = 1.0 + 2(c2y)
    yB_A[3] = 2(c2z + qs*qx)

    return nothing
end

"""
    q_toz!(zB_A, q_AB)

Compute the z axis of frame B projected in frame A.
"""
@inline function q_toz!(zB_A, q_AB)
    qs, qx, qy, qz = q_AB

    c2x = qz*qx
    c2y = qz*qy
    c2z = -qx*qx - qy*qy

    zB_A[1] = 2(c2x + qs*qy)
    zB_A[2] = 2(c2y - qs*qx)
    zB_A[3] = 1.0 + 2(c2z)

    return nothing
end

# This function determines the convention used by a quaternion multiplication function in
# terms of: 1) Scalar-first or scalar-last, 2) Right-handed or left-handed algebra
# The multiplication function shall provide q*p = fmult(q, p)
function q_testConvention(fmult)
    q = [1.0; 2.0; 3.0; 4.0]
    p = [-4.0; 2.0; -3.0; 1.0]
    q_testConvention(fmult(q, p), q, p)
end

function q_testConvention(qp, q, p)

    # out = q*p
    function q_multiplyUniversal(q, p; scalarFirst=true, right=true)
        qq = copy(q);
        pp = copy(p)
        if !scalarFirst
            qq = [q[4]; q[1:3]]
            pp = [p[4]; p[1:3]]
        end
        qqpp = right ? qp_right(qq, pp) : qp_left(qq, pp)
        if scalarFirst
            ;
            return qqpp;
        end
        return [qqpp[2:4]; qqpp[1]]
    end

    function qp_right(q, p)
        qs, qx, qy, qz = q
        return [
            +qs -qx -qy -qz;
            +qx +qs -qz +qy;
            +qy +qz +qs -qx;
            +qz -qy +qx +qs
        ]*p
    end

    function qp_left(q, p)
        qs, qx, qy, qz = q
        return [
            +qs -qx -qy -qz;
            +qx +qs +qz -qy;
            +qy -qz +qs +qx;
            +qz +qy -qx +qs
        ]*p
    end

    qerr(a, b) = min(norm(a - b), norm(a + b))
    qp1 = q_multiplyUniversal(q, p; scalarFirst=true, right=true)
    qp2 = q_multiplyUniversal(q, p; scalarFirst=false, right=true)
    qp3 = q_multiplyUniversal(q, p; scalarFirst=true, right=false)
    qp4 = q_multiplyUniversal(q, p; scalarFirst=false, right=false)
    errv = [qerr(qp1, qp), qerr(qp2, qp), qerr(qp3, qp), qerr(qp4, qp)]

    msgs = [
        "Scalar first, right-handed (ij = k)",
        "Scalar last, right-handed (ij = k)",
        "Scalar first, left-handed (ij = -k)",
        "Scalar-last, left-handed (ij = -k)",
    ]
    println("Quaternion convention: \n"*msgs[findmin(errv)[2]])
end
