"""
    q_AC = q_multiply(q_AB, q_BC)

Mutliply the two input quaternions as follows:
```math
q_{AC} = q_{AB} ⊗ q_{BC}
```
"""
@inline q_multiply(q_AB, q_BC) = collect(q_multiplyCore(q_AB, q_BC))
@inline q_multiply(q_AB::SVector{4}, q_BC::SVector{4}) = SVector(q_multiplyCore(q_AB, q_BC))
@inline q_multiply(q_AB::Quaternion, q_BC::Quaternion) = q_multiplyCore(q_AB, q_BC)
@inline q_multiply(q_AB::Quaternion, q_BC::AbstractVector) = q_multiplyCore(q_AB, q_BC)
@inline q_multiply(q_AB::AbstractVector, q_BC::Quaternion) = q_multiplyCore(q_AB, q_BC)
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
    qOut = copy(q[1])
    for i in 2:lastindex(q)
        q_multiply!(qOut, q[i])
    end
    return qOut
end

function q_multiplyn!(qOut, q...)
    qOut .= q[1]
    for i in 2:lastindex(q)
        q_multiply!(qOut, q[i])
    end
    return qOut
end

"""
    R_AB = q_toDcm(q_AB)

Translate the input unitary quaternion into a transformation matrix.
```math
R_{AB}(q_{AB}) = I + 2qₛ[qᵥ×] + 2[qᵥ×]²
```
"""
@inline function q_toDcm(q_AB)
    R_AB = Matrix{eltype(q_AB)}(undef, 3, 3)
    return q_toDcm!(R_AB, q_AB)
end

@inline q_toDcm(q_AB::SVector{4}) = SMatrix{3, 3}(q_toDcmCore(q_AB))

@inline function q_toDcm(q_AB::Quaternion)
    r11, r21, r31, r12, r22, r32, r13, r23, r33 = q_toDcmCore(q_AB)
    return [r11 r12 r13; r21 r22 r23; r31 r32 r33]
end

@inline function q_toDcm!(R_AB, q_AB)     # R_AB from q_AB
    R_AB[1, 1], R_AB[2, 1], R_AB[3, 1],
    R_AB[1, 2], R_AB[2, 2], R_AB[3, 2],
    R_AB[1, 3], R_AB[2, 3], R_AB[3, 3] = q_toDcmCore(q_AB)
    return R_AB
end

"""
    q_AB = q_fromDcm(R_AB)

Translate the input rotation matrix into a unitary quaternion.
"""
@inline q_fromDcm(R_AB) = collect(q_fromDcmCore(R_AB))
@inline q_fromDcm!(q_AB, R_AB) = q_AB .= q_fromDcmCore(R_AB)
@inline q_fromDcm(R_AB::SMatrix{3, 3}) = SVector{4}(q_fromDcmCore(R_AB))

"""
    q_AB = q_fromAxes(xB_A, yB_A, zB_A)

Compute the attitude quaternion given as input the axes of a reference frame.
"""
# @inline q_fromAxes(xB_A, yB_A, zB_A) = q_fromDcm(dcm_fromAxes(xB_A, yB_A, zB_A))

@inline function q_fromAxes(xB_A, yB_A, zB_A)
    q_AB = Vector{eltype(yB_A)}(undef, 4)
    return q_fromAxes!(q_AB, xB_A, yB_A, zB_A)
end

@inline function q_fromAxes!(q_AB, xB_A, yB_A, zB_A)
    r11, r21, r31 = xB_A
    r12, r22, r32 = yB_A
    r13, r23, r33 = zB_A
    R_AB = (r11, r21, r31, r12, r22, r32, r13, r23, r33)
    return q_AB .= q_fromDcmCore(R_AB)
end

@inline function q_fromAxes(xB_A::Vector3, yB_A::Vector3, zB_A::Vector3)
    r11, r12, r13 = xB_A
    r21, r22, r23 = yB_A
    r31, r32, r33 = zB_A
    R_AB = (r11, r21, r31, r12, r22, r32, r13, r23, r33)
    return q_fromDcmCore(R_AB)
end

"""
    q_random()

Generate a random unitary quaternion.
"""
@inline q_random() = normalize!(randn(4))
@inline q_random(rng) = normalize!(randn(rng), 4)

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
@inline q_transformVector(q_AB, v_B) = collect(q_transformVectorCore(q_AB, v_B))
@inline q_transformVector(q_AB::SVector{4}, v_B::SVector{3}) = SVector{3}(q_transformVectorCore(q_AB, v_B))
@inline q_transformVector(q_AB::AbstractVector, v_B::Vector3) = q_transformVectorCore(q_AB, v_B)
@inline q_transformVector(q_AB::Quaternion, v_B::AbstractVector) = q_transformVectorCore(q_AB, v_B)
@inline q_transformVector(q_AB::Quaternion, v_B::Vector3) = q_transformVectorCore(q_AB, v_B)
@inline q_transformVector!(v_A, q_AB, v_B) = v_A .= q_transformVectorCore(q_AB, v_B)

"""
    q_transformVectorT!(v_A, q_BA, v_B)

Project the vector v from frame B into frame A, using q_BA.
"""
@inline q_transformVectorT(q_BA, v_B) = collect(q_transformVectorCore(q_transposeCore(q_BA), v_B))
@inline q_transformVectorT(q_BA::AbstractVector, v_B::Vector3) = q_transformVectorCore(q_transposeCore(q_BA), v_B)
@inline q_transformVectorT(q_BA::Quaternion, v_B::AbstractVector) = q_transformVectorCore(q_transposeCore(q_BA), v_B)
@inline q_transformVectorT(q_BA::Quaternion, v_B::Vector3) = q_transformVectorCore(q_transposeCore(q_BA), v_B)
@inline q_transformVectorT(q_BA::SVector{4}, v_B::SVector{3}) = SVector{3}(q_transformVectorCore(q_transposeCore(q_BA), v_B))
@inline q_transformVectorT!(v_A, q_BA, v_B) = v_A .= q_transformVectorCore(q_transposeCore(q_BA), v_B)

"""
    q' = q_transpose(q)

Transpose the input quaternion.
"""
@inline q_transpose(q) = collect(q_transposeCore(q))
@inline q_transpose(q::SVector{4}) = SVector{4}(q_transposeCore(q))
@inline q_transpose(q::Quaternion) = q_transposeCore(q)
@inline q_transpose!(q) = q .= q_transposeCore(q)
@inline q_transpose!(qt, q) = qt .= q_transposeCore(q)

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
@inline q_derivative(q_AB, ωAB_B) = collect(q_derivativeCore(q_AB, ωAB_B))
@inline q_derivative(q_AB::SVector{4}, ωAB_B::SVector{3}) = SVector{4}(q_derivativeCore(q_AB, ωAB_B))
@inline q_derivative(q_AB::Quaternion, ωAB_B) = q_derivativeCore(q_AB, ωAB_B)
@inline q_derivative(q_AB, ωAB_B::Vector3) = q_derivativeCore(q_AB, ωAB_B)
@inline q_derivative!(dq_AB, q_AB, ωAB_B) = dq_AB .= q_derivativeCore(q_AB, ωAB_B)

@inline function q_rotx(θ)
    sθ, cθ = sincos(θ / 2)
    return [cθ; sθ; 0; 0]
end

@inline function q_rotx!(q, θ)
    q[2], q[1] = sincos(θ / 2)
    q[3] = q[4] = 0
    return q
end

@inline function q_roty(θ)
    sθ, cθ = sincos(θ / 2)
    return [cθ; 0; sθ; 0]
end

@inline function q_roty!(q, θ)
    q[3], q[1] = sincos(θ / 2)
    q[2] = q[4] = 0
    return q
end

@inline function q_rotz(θ)
    sθ, cθ = sincos(θ / 2)
    return [cθ; 0; 0; sθ]
end

@inline function q_rotz!(q, θ)
    q[4], q[1] = sincos(θ / 2)
    q[2] = q[3] = 0
    return q
end

"""
    q_AB = q_fromAxisAngle(u, θ_AB)

Compute the unitary quaternion given as input an axis-angle representation.
"""
q_fromAxisAngle(u, θ) = iszero(u) ? q_identity() : [cos(0.5θ); sin(0.5θ)*normalize(u)]

@inline function q_fromAxisAngle(idx::Int, θ)
    sθ, cθ = sincos(θ / 2)
    q = [cθ; 0; 0; 0]
    q[idx + 1] = sθ
    return q
end

function q_fromAxisAngle!(q, idx::Int, θ)
    sθ, cθ = sincos(θ / 2)
    q .= 0
    q[1] = cθ
    q[idx + 1] = sθ
    return q
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

@inline function q_toAxes(q_AB)
    R = q_toDcm(q_AB)
    return R[1, :], R[2, :], R[3, :]    # xB_A, yB_A, zB_A
end


"""
    q⁻¹ = q_inverse(q)

Compute the inverse of the input quaternion.
"""
@inline q_inverse(q) = collect(q_inverseCore(q))
@inline q_inverse(q::Vector3) = q_inverseCore(q)
@inline q_inverse(q::SVector{4}) = SVector{4}(q_inverseCore(q))

@inline q_toRv(q) = collect(q_toRvCore(q))
@inline q_toRv(q::Vector3)= q_toRvCore(q)
@inline q_toRv(q::SVector{4}) = SVector{4}(q_toRvCore(q))

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
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:xyz}) = q_fromSequence((θ1, θ2, θ3), (1, 2, 3))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:xzy}) = q_fromSequence((θ1, θ2, θ3), (1, 3, 2))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:yxz}) = q_fromSequence((θ1, θ2, θ3), (2, 1, 3))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:yzx}) = q_fromSequence((θ1, θ2, θ3), (2, 3, 1))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:zxy}) = q_fromSequence((θ1, θ2, θ3), (3, 1, 2))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:zyx}) = q_fromSequence((θ1, θ2, θ3), (3, 2, 1))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:xyx}) = q_fromSequence((θ1, θ2, θ3), (1, 2, 1))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:xzx}) = q_fromSequence((θ1, θ2, θ3), (1, 3, 1))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:yxy}) = q_fromSequence((θ1, θ2, θ3), (2, 1, 2))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:yzy}) = q_fromSequence((θ1, θ2, θ3), (2, 3, 2))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:zxz}) = q_fromSequence((θ1, θ2, θ3), (3, 1, 3))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:zyz}) = q_fromSequence((θ1, θ2, θ3), (3, 2, 3))

@inline function q_fromSequence(θ, sequence=(3, 2, 1))
    q = q_fromAxisAngle(sequence[1], θ[1])
    for i in 2:lastindex(θ)
        q_multiply!(q, q_fromAxisAngle(sequence[i], θ[i]))
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
    return [k*cos(qvn); qx*qvnn; qy*qvnn; qz*qvnn]
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
    θerr = Vector{eltype(q)}(undef, 3)
    return q_attitudeError!(θerr, qNominal, q)
end

@inline function q_attitudeError!(θerr, qNominal, q)
    ps, px, py, pz = q
    qs, qx, qy, qz = qNominal

    imax = findmax(abs, qNominal)[2]
    sgn = sign(qNominal[imax]) == sign(q[imax]) ? 2.0 : -2.0

    θerr[1] = (-px*qs + ps*qx + pz*qy - py*qz)*sgn
    θerr[2] = (-py*qs - pz*qx + ps*qy + px*qz)*sgn
    θerr[3] = (-pz*qs + py*qx - px*qy + ps*qz)*sgn

    return θerr  # 2q_multiply(q_transpose(sgn*q), qNominal)[2:4]
end

"""
    xB_A = q_tox(q_AB)

Compute the x axis of frame B projected in frame A.
"""
@inline function q_tox(q_AB)
    out = Vector{eltype(q_AB)}(undef, 3)
    return q_tox!(out, q_AB)
end

"""
    yB_A = q_toy(q_AB)

Compute the y axis of frame B projected in frame A.
"""
@inline function q_toy(q_AB)
    out = Vector{eltype(q_AB)}(undef, 3)
    return q_toy!(out, q_AB)
end

"""
    zB_A = q_toz(q_AB)

Compute the z axis of frame B projected in frame A.
"""
@inline function q_toz(q_AB)
    out = Vector{eltype(q_AB)}(undef, 3)
    return q_toz!(out, q_AB)
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

    return xB_A
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

    return yB_A
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

    return zB_A
end

function q_fromZaxis(zB_A)
    zB_A = normalize(zB_A)
    return q_fromDcm([nullspace(zB_A') zB_A])
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
        qq = copy(q)
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
