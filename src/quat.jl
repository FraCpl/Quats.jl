@inline Quats.Quat(q::AbstractVector) = q[1], q[2], q[3], q[4]
@inline q_zero() = 1.0, 0.0, 0.0, 0.0
@inline Base.:*(s::Number, a::Quats.Quat) = s*a[1], s*a[2], s*a[3], s*a[4]              # Scalar multiplication
@inline Base.:*(a::Quats.Quat, s::Number) = s*a[1], s*a[2], s*a[3], s*a[4]              # Scalar multiplication
@inline Base.:/(a::Quats.Quat, s::Number) = a[1]/s, a[2]/s, a[3]/s, a[4]/s              # Scalar division
@inline Base.:*(a::Quats.Quat, b::Quats.Quat) = q_multiplyCore(a, b)                    # Quaternion multiplication
@inline Base.:*(a::Quats.Quat, b::Quats.Vec3) = q_transformVector(a, b)                 # Vector transformation
@inline LinearAlgebra.norm(q::Quats.Quat) = hypot(q[1], q[2], q[3], q[4])
@inline LinearAlgebra.normalize(q::Quats.Quat) = q / norm(q)
@inline LinearAlgebra.conj(q::Quats.Quat) = q[1], -q[2], -q[3], -q[4]
@inline LinearAlgebra.adjoint(q::Quats.Quat) = q[1], -q[2], -q[3], -q[4]

@inline toVector(q::Quats.Quat) = [q[1], q[2], q[3], q[4]]
@inline Base.copyto!(v::AbstractVector, q::Quats.Quat) = v .= q

@inline randn4() = randn(), randn(), randn(), randn()
@inline randn4(rng) = randn(rng), randn(rng), randn(rng), randn(rng)
@inline rand4() = rand(), rand(), rand(), rand()
@inline rand4(rng) = rand(rng), rand(rng), rand(rng), rand(rng)

# QUATERNION MULTIPLICATION -------------------------------------------------------------- #
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

"""
    q_AC = q_multiply(q_AB, q_BC)

Mutliply the two input quaternions as follows:
```math
q_{AC} = q_{AB} ⊗ q_{BC}
```
"""
@inline q_multiply(q_AB::AbstractVector, q_BC::AbstractVector) = toVector(q_multiplyCore(q_AB, q_BC))
@inline q_multiply(q_AB, q_BC) = q_multiplyCore(q_AB, q_BC)
@inline q_multiply!(q_AC, q_AB, q_BC) = q_AC .= q_multiplyCore(q_AB, q_BC)
@inline q_multiply!(q1, q2) = q1 .= q_multiplyCore(q1, q2)      # q1 is modified in-place

# Auiliary multiplication functions: these are included to allow high-speed in-place
# multiplication operations without the need of allocating/computing a quaternion transpose
@inline q_multiplyT1!(q_AC, q_BA, q_BC) = q_AC .= q_multiplyCore(q_transposeCore(q_BA), q_BC)
@inline q_multiplyT2!(q_AC, q_AB, q_CB) = q_AC .= q_multiplyCore(q_AB, q_transposeCore(q_CB))
@inline q_multiplyT12!(q_AC, q_BA, q_CB) = q_AC .= q_transposeCore(q_multiplyCore(q_CB, q_BA))

"""
    q = q_multiplyn(q1, q2, q3, ...)

```math
    q = q₁ ⊗ q₂ ⊗ q₃ ⊗ ...
```
"""
@inline function q_multiplyn(q...)
    s, x, y, z = q[1]
    qOut = (s, x, y, z)
    for i in 2:lastindex(q)
        qOut = q_multiplyCore(qOut, q[i])
    end
    return qOut
end

q_multiplyn!(qOut, q...) = qOut .= q_multiplyn(q...)

# QUATERNION TO DCM ---------------------------------------------------------------------- #
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

"""
    R_AB = q_toDcm(q_AB)

Translate the input unitary quaternion into a transformation matrix.
```math
R_{AB}(q_{AB}) = I + 2qₛ[qᵥ×] + 2[qᵥ×]²
```
"""
@inline q_toDcm(q_AB::AbstractVector) = toMatrix(q_toDcmCore(q_AB))
@inline q_toDcm(q_AB) = q_toDcmCore(q_AB)
@inline q_toDcm!(R_AB, q_AB) = toMatrix!(R_AB, q_toDcmCore(q_AB))

# QUATERNION FROM DCM -------------------------------------------------------------------- #
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

"""
    q_AB = q_fromDcm(R_AB)

Translate the input rotation matrix into a unitary quaternion.
"""
@inline q_fromDcm(R_AB::AbstractMatrix) = toVector(q_fromDcmCore(R_AB))
@inline q_fromDcm(R_AB) = q_fromDcmCore(R_AB)
@inline q_fromDcm!(q_AB, R_AB) = q_AB .= q_fromDcmCore(R_AB)

"""
    q_random()

Generate a random unitary quaternion.
"""
@inline q_random() = normalize(randn4())
@inline q_random(rng) = normalize(randn4(rng))

# VECTOR TRANSFORMATION ------------------------------------------------------------------ #
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

"""
    v_A = q_transformVector(q_AB, v_B)

Project the vector v from frame B into frame A.
"""
@inline q_transformVector(q_AB::AbstractVector, v_B::AbstractVector) = toVector(q_transformVectorCore(q_AB, v_B))
@inline q_transformVector(q_AB, v_B) = q_transformVectorCore(q_AB, v_B)
@inline q_transformVector!(v_A, q_AB, v_B) = v_A .= q_transformVectorCore(q_AB, v_B)

"""
    q_transformVectorT!(v_A, q_BA, v_B)

Project the vector v from frame B into frame A, using q_BA.
"""
@inline q_transformVectorT(q_BA::AbstractVector, v_B::AbstractVector) = toVector(q_transformVectorCore(q_transposeCore(q_BA), v_B))
@inline q_transformVectorT(q_BA, v_B) = q_transformVectorCore(q_transposeCore(q_BA), v_B)
@inline q_transformVectorT!(v_A, q_BA, v_B) = v_A .= q_transformVectorCore(q_transposeCore(q_BA), v_B)

# QUATERNION TRANSPOSE ------------------------------------------------------------------- #
@inline q_transposeCore(q) = q[1], -q[2], -q[3], -q[4]

"""
    q' = q_transpose(q)

Transpose the input quaternion.
"""
@inline q_transpose(q::AbstractVector) = toVector(q_transposeCore(q))
@inline q_transpose(q) = q_transposeCore(q)
@inline q_transpose!(qt, q) = qt .= q_transposeCore(q)

# QUATERNION DERIVATIVE ------------------------------------------------------------------ #
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
@inline q_derivative(q_AB::AbstractVector, ωAB_B::AbstractVector) = toVector(q_derivativeCore(q_AB, ωAB_B))
@inline q_derivative(q_AB, ωAB_B) = q_derivativeCore(q_AB, ωAB_B)
@inline q_derivative!(dq_AB, q_AB, ωAB_B) = dq_AB .= q_derivativeCore(q_AB, ωAB_B)

# SIMPLE ROTATIONS ----------------------------------------------------------------------- #
@inline function q_rotx(θ::T) where {T}
    sθ, cθ = sincos(θ / 2)
    o = zero(T)
    return cθ, sθ, o, o
end
@inline q_rotx!(q, θ) = q .= q_rotx(θ)

@inline function q_roty(θ::T) where {T}
    sθ, cθ = sincos(θ / 2)
    o = zero(T)
    return cθ, o, sθ, o
end
@inline q_roty!(q, θ) = q .= q_roty(θ)

@inline function q_rotz(θ::T) where {T}
    sθ, cθ = sincos(θ / 2)
    o = zero(T)
    return cθ, o, o, sθ
end
@inline q_rotz!(q, θ) = q .= q_rotz(θ)

"""
    q_AB = q_fromAxisAngle(u, θ_AB)

Compute the unitary quaternion given as input an axis-angle representation.
"""
q_fromAxisAngle(u, θ) = iszero(u) ? q_zero() : [cos(0.5θ); sin(0.5θ)*normalize(u)]

function q_fromAxisAngle(u::Vec3, θ)
    iszero(u) && return (1.0, 0.0, 0.0, 0.0)
    x, y, z = normalize(u)
    sθ, cθ = sincos(θ / 2)
    return cθ, x*sθ, y*sθ, z*sθ
end

function q_fromAxisAngle!(q, u, θ)
    uNormSq = 0.0
    @inbounds for i in 1:3
        uNormSq += u[i] * u[i]
    end

    if uNormSq == 0
        q .= 0.0
        q[1] = 1.0
        return q
    end

    st, q[1] = sincos(θ / 2)
    st = st / sqrt(uNormSq)
    @inbounds for i in 1:3
        q[i + 1] = st * u[i]
    end
    return q
end

@inline function q_toAxisAngle(q)
    qs, qx, qy, qz = q
    nqv = sqrt(qx*qx + qy*qy + qz*qz)
    return [qx/nqv; qy/nqv; qz/nqv], 2atan(nqv, qs)
end

@inline function q_fromAxesCore(xB_A, yB_A, zB_A)
    r11, r21, r31 = xB_A
    r12, r22, r32 = yB_A
    r13, r23, r33 = zB_A
    return q_fromDcmCore((r11, r21, r31, r12, r22, r32, r13, r23, r33)) # q_AB
end

"""
    q_AB = q_fromAxes(xB_A, yB_A, zB_A)

Compute the attitude quaternion given as input the axes of a reference frame.
"""
@inline q_fromAxes(xB_A::AbstractVector, yB_A::AbstractVector, zB_A::AbstractVector) = toVector(q_fromAxesCore(xB_A, yB_A, zB_A))
@inline q_fromAxes(xB_A, yB_A, zB_A) = q_fromAxesCore(xB_A, yB_A, zB_A)
@inline q_fromAxes!(q_AB, xB_A, yB_A, zB_A) = q_AB .= q_fromAxesCore(xB_A, yB_A, zB_A)

@inline function q_toAxes(q_AB)
    r11, r21, r31, r12, r22, r32, r13, r23, r33 = q_toDcmCore(q_AB)
    return (r11, r21, r31), (r12, r22, r32), (r13, r23, r33)    # xB_A, yB_A, zB_A
end

# QUATERNION INVERSE --------------------------------------------------------------------- #
function q_inverseCore(q)
    qs, qx, qy, qz = q
    nq2 = qs * qs + qx * qx + qy * qy + qz * qz
    return qs / nq2, -qx / nq2, -qy / nq2, -qz / nq2
end

"""
    q⁻¹ = q_inverse(q)

Compute the inverse of the input quaternion.
"""
@inline q_inverse(q::AbstractVector) = toVector(q_inverseCore(q))
@inline q_inverse(q) = q_inverseCore(q)


@inline function q_toRvCore(q)
    qs, qx, qy, qz = q
    nqv = sqrt(qx*qx + qy*qy + qz*qz)
    f = 2*atan(nqv, qs) / ifelse(iszero(nqv), one(nqv), nqv)
    return f*qx, f*qy, f*qz
end

@inline q_toRv(q::AbstractVector) = toVector(q_toRvCore(q))
@inline q_toRv(q)= q_toRvCore(q)

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
function q_rate(Δt, q1_AB, q2_AB)
    dq_AB = q_multiply(q_transpose(q1_AB), q2_AB)
    return q_toRv(dq_AB) ./ Δt
end

# https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
# CAUTION: Only :zyx and :xyz sequences are implemented for now
# θz, θy, θx = q_toEuler(q, :zyx)
# θx, θz, θy = q_toEuler(q, :xzy)
@inline q_toEuler(q, s::Symbol=:zyx) = q_toEuler(q, Val(s))

@inline function q_toEuler(q, ::Val{:zyx})
    qs, qx, qy, qz = q

    t1 = 2 * (qs*qx + qy*qz)
    t2 = 1 - 2 * (qx*qx + qy*qy)
    θx = atan(t1, t2)

    t3 = 2 * (qs*qy - qz*qx)
    t3 = clamp(t3, -1.0, 1.0)
    θy = asin(t3)

    t4 = 2 * (qs*qz + qx*qy)
    t5 = 1 - 2 * (qy*qy + qz*qz)
    θz = atan(t4, t5)

    return θz, θy, θx
end

function q_toEuler(q, ::Val{:xyz})
    qs, qx, qy, qz = q

    t2 = 2 * (qx*qz + qy*qs)
    t2 = clamp(t2, -1.0, 1.0)
    θy = asin(t2)

    if abs(t2) ≥ 1 - 10*eps(eltype(q))
        θx = zero(t2)
        θz = sign(t2)*2*atan(qx, qs)
    else
        θx = atan(2 * (qx*qs - qy*qz), qs^2 - qx^2 - qy^2 + qz^2)
        θz = atan(2 * (qz*qs - qx*qy), qs^2 + qx^2 - qy^2 - qz^2)
    end

    return θx, θy, θz
end

@inline function q_toEuler(q, ::Val{:xzy})
    qs, qx, qy, qz = q

    t2 = 2 * (qs * qz - qx * qy)
    t2 = clamp(t2, -1.0, 1.0)
    θz = asin(t2)

    if abs(t2) ≥ 1 - 10*eps(eltype(q))
        θx = zero(t2)
        θy = sign(t2)*2*atan(qx, qs)
    else
        θx = atan(2 * (qx*qs + qy*qz), 1 - 2*(qx^2 + qz^2))
        θy = atan(2 * (qy*qs + qx*qz), 1 - 2*(qy^2 + qz^2))
    end

    return θx, θz, θy
end

@inline function q_toEuler(q, ::Val{:yxz})
    qs, qx, qy, qz = q

    t2 = 2*(qs*qx - qy*qz)
    t2 = clamp(t2, -1.0, 1.0)
    θx = asin(t2)

    if abs(t2) ≥ 1 - 10*eps(eltype(q))
        θy = zero(t2)
        θz = sign(t2)*2*atan(qy, qs)
    else
        θy = atan(2*(qy*qs + qx*qz), 1 - 2*(qx^2 + qy^2))
        θz = atan(2*(qz*qs + qx*qy), 1 - 2*(qx^2 + qz^2))
    end

    return θy, θx, θz
end

# @inline function q_toEuler(q, ::Val{:yzx})
#     qs, qx, qy, qz = q

#     t2 = 2*(qy*qx + qs*qz)
#     t2 = clamp(t2, -1.0, 1.0)
#     θz = asin(t2)

#     if abs(t2) ≥ 1 - 10*eps(eltype(q))
#         θy = zero(t2)
#         θx = sign(t2)*2*atan(qy, qs)
#     else
#         θy = atan(-2*(qx*qz - qs*qy), 1 - 2*(qy^2 + qz^2))
#         θx = atan(-2*(qx*qy - qs*qx), 1 - 2*(qx^2 + qz^2))
#     end

#     return θy, θz, θx
# end

# (θ1=θ_AB, θ2=θ_BC, θ3=θ_CD) --> q_AD
@inline q_fromEuler(θ1, θ2, θ3, s::Symbol=:zyx) = q_fromEuler(θ1, θ2, θ3, Val(s))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:xyz}) = q_multiplyn(q_rotx(θ1), q_roty(θ2), q_rotz(θ3))#q_fromSequence((θ1, θ2, θ3), (1, 2, 3))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:xzy}) = q_multiplyn(q_rotx(θ1), q_rotz(θ2), q_roty(θ3))#q_fromSequence((θ1, θ2, θ3), (1, 3, 2))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:yxz}) = q_multiplyn(q_roty(θ1), q_rotx(θ2), q_rotz(θ3))#q_fromSequence((θ1, θ2, θ3), (2, 1, 3))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:yzx}) = q_multiplyn(q_roty(θ1), q_rotz(θ2), q_rotx(θ3))#q_fromSequence((θ1, θ2, θ3), (2, 3, 1))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:zxy}) = q_multiplyn(q_rotz(θ1), q_rotx(θ2), q_roty(θ3))#q_fromSequence((θ1, θ2, θ3), (3, 1, 2))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:zyx}) = q_multiplyn(q_rotz(θ1), q_roty(θ2), q_rotx(θ3))#q_fromSequence((θ1, θ2, θ3), (3, 2, 1))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:xyx}) = q_multiplyn(q_rotx(θ1), q_roty(θ2), q_rotx(θ3))#q_fromSequence((θ1, θ2, θ3), (1, 2, 1))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:xzx}) = q_multiplyn(q_rotx(θ1), q_rotz(θ2), q_rotx(θ3))#q_fromSequence((θ1, θ2, θ3), (1, 3, 1))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:yxy}) = q_multiplyn(q_roty(θ1), q_rotx(θ2), q_roty(θ3))#q_fromSequence((θ1, θ2, θ3), (2, 1, 2))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:yzy}) = q_multiplyn(q_roty(θ1), q_rotz(θ2), q_roty(θ3))#q_fromSequence((θ1, θ2, θ3), (2, 3, 2))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:zxz}) = q_multiplyn(q_rotz(θ1), q_rotx(θ2), q_rotz(θ3))#q_fromSequence((θ1, θ2, θ3), (3, 1, 3))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:zyz}) = q_multiplyn(q_rotz(θ1), q_roty(θ2), q_rotz(θ3))#q_fromSequence((θ1, θ2, θ3), (3, 2, 3))


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
    return (k*cos(qvn), qx*qvnn, qy*qvnn, qz*qvnn)
end

"""
    ql = q_log(q)

Compute the logarithm of the input quaternion.
"""
@inline function q_log(q, n=1)
    qs, qx, qy, qz = q
    vNorm = sqrt(qx*qx + qy*qy + qz*qz)
    if iszero(vNorm)
        vNorm = 1.0
    end
    nq = norm(q)
    vNorm = acos(qs / nq) / vNorm
    return (n*log(nq), n*qx*vNorm, n*qy*vNorm, n*qz*vNorm)
end

@inline q_power(q, n) = q_exp(q_log(q, n))

"""
    q = q_slerp(q0, q1, τ)

Compute the spherical linear interpolation between two quaternions at τ, where q0 = q[τ=0], q1 = q[τ=1], and q = q[τ].
"""
@inline q_slerp(q0, q1, τ) = q_multiply(q0, q_power(q_multiply(q_transpose(q0), q1), τ))

@inline function q_interp(t, q, ti)
    id0 = findlast(t .≤ ti)
    id1 = min(id0+1, length(t))
    return q_slerp(q[id0], q[id1], (ti - t[id0])/(t[id1] - t[id0]))
end


@inline function q_attitudeErrorCore(qNominal, q)
    ps, px, py, pz = q
    qs, qx, qy, qz = qNominal

    imax = findmax(abs, qNominal)[2]
    sgn = copysign(2.0, qNominal[imax] * q[imax])
    # sgn = sign(qNominal[imax]) == sign(q[imax]) ? 2 : -2

    θx = (-px*qs + ps*qx + pz*qy - py*qz)*sgn
    θy = (-py*qs - pz*qx + ps*qy + px*qz)*sgn
    θz = (-pz*qs + py*qx - px*qy + ps*qz)*sgn

    return θx, θy, θz  # 2q_multiply(q_transpose(sgn*q), qNominal)[2:4]
end


@inline q_attitudeError!(θerr, qNominal, q) = θerr .= q_attitudeErrorCore(qNominal, q)
@inline q_attitudeError(qNominal, q) = q_attitudeErrorCore(qNominal, q)

function q_fromZaxis(zB_A)
    zB_A = normalize(zB_A)
    _, i = findmin(abs, zB_A)
    v = i == 1 ? Quats.ex : i == 2 ? Quats.ey : Quats.ez
    xB_A = normalize(cross(zB_A, v))
    yB_A = cross(zB_A, xB_A)
    return q_fromAxes(xB_A, yB_A, zB_A)     # q_AB
end

"""
    xB_A = q_tox(q_AB)

Compute the x axis of frame B projected in frame A.
"""
@inline function q_tox(q_AB)
    qs, qx, qy, qz = q_AB

    c2x = -qy*qy - qz*qz
    c2y = qx*qy
    c2z = qx*qz

    # Compute xB_A
    x = 1.0 + 2*c2x
    y = 2(c2y + qs*qz)
    z = 2(c2z - qs*qy)
    return x, y, z
end

"""
    yB_A = q_toy(q_AB)

Compute the y axis of frame B projected in frame A.
"""
@inline function q_toy(q_AB)
    qs, qx, qy, qz = q_AB

    c2x = qy*qx
    c2y = -qz*qz - qx*qx
    c2z = qy*qz

    # Compute yB_A
    x = 2(c2x - qs*qz)
    y = 1.0 + 2(c2y)
    z = 2(c2z + qs*qx)
    return x, y, z
end

"""
    zB_A = q_toz(q_AB)

Compute the z axis of frame B projected in frame A.
"""
@inline function q_toz(q_AB)
    qs, qx, qy, qz = q_AB

    c2x = qz*qx
    c2y = qz*qy
    c2z = -qx*qx - qy*qy

    # Copmute zB_A
    x = 2(c2x + qs*qy)
    y = 2(c2y - qs*qx)
    z = 1 + 2(c2z)
    return x, y, z
end

# DISPLAY FUNCTIONS ---------------------------------------------------------------------- #
function Base.show(io::IO, ::MIME"text/plain", q::Quat{T}) where T
    println(io, "Quats.Quat{$T}:")
    for i in 1:4
        print(io, " ")
        x = q[i]
        !signbit(x) && print(io, " ")
        show(io, round(x, digits=16))
        i < 4 && println(io)
    end
end
