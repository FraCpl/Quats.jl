"""
    q_identity

Get identity quaternion, with scalar component equal to 1 and
vector components equal to zero.
"""
q_identity() = [1.0, 0.0, 0.0, 0.0]

"""
    q_AC = q_multiply(q_AB, q_BC)

Mutliply the two input quaternions as follows:
```math
q_{AC} = q_{AB} ⊗ q_{BC}
```
"""
@inline q_multiply(q_AB, q_BC) = collect(q_multiplyCore(q_AB, q_BC))
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
@inline q_toDcm(q_AB::AbstractVector) = toMatrix(q_toDcmCore(q_AB))
@inline q_toDcm(q_AB::Quaternion) = q_toDcmCore(q_AB)
@inline q_toDcm!(R_AB::AbstractMatrix, q_AB) = toMatrix!(R_AB, q_toDcmCore(q_AB))

"""
    q_AB = q_fromDcm(R_AB)

Translate the input rotation matrix into a unitary quaternion.
"""
@inline q_fromDcm(R_AB::AbstractMatrix) = collect(q_fromDcmCore(R_AB))
@inline q_fromDcm(R_AB::Matrix3x3) = q_fromDcmCore(R_AB)
@inline q_fromDcm!(q_AB::AbstractVector, R_AB) = q_AB .= q_fromDcmCore(R_AB)

"""
    q_AB = q_fromAxes(xB_A, yB_A, zB_A)

Compute the attitude quaternion given as input the axes of a reference frame.
"""
@inline q_fromAxes!(q_AB::AbstractVector, xB_A, yB_A, zB_A) = q_AB .= q_fromAxesCore(xB_A, yB_A, zB_A)
@inline q_fromAxes(xB_A::Vector3, yB_A::Vector3, zB_A::Vector3) = q_fromAxesCore(xB_A, yB_A, zB_A)
@inline q_fromAxes(xB_A::AbstractVector, yB_A::AbstractVector, zB_A::AbstractVector) = toVector(q_fromAxesCore(xB_A, yB_A, zB_A))

"""
    q_random()

Generate a random unitary quaternion.
"""
@inline q_random() = normalize!(randn(4))
@inline q_random(rng) = normalize!(randn(rng, 4))

"""
    v_A = q_transformVector(q_AB, v_B)

Project the vector v from frame B into frame A.
"""
@inline q_transformVector(q_AB, v_B) = collect(q_transformVectorCore(q_AB, v_B))
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
@inline q_transformVectorT!(v_A, q_BA, v_B) = v_A .= q_transformVectorCore(q_transposeCore(q_BA), v_B)

"""
    q' = q_transpose(q)

Transpose the input quaternion.
"""
@inline q_transpose(q::AbstractVector) = collect(q_transposeCore(q))
@inline q_transpose(q::Quaternion) = q_transposeCore(q)
@inline q_transpose!(q::AbstractVector) = q .= q_transposeCore(q)
@inline q_transpose!(qt::AbstractVector, q::Quaternion) = qt .= q_transposeCore(q)
@inline q_transpose!(qt::AbstractVector, q::AbstractVector) = qt .= q_transposeCore(q)

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
@inline q_derivative(q_AB::Quaternion, ωAB_B) = q_derivativeCore(q_AB, ωAB_B)
@inline q_derivative(q_AB, ωAB_B::Vector3) = q_derivativeCore(q_AB, ωAB_B)
@inline q_derivative!(dq_AB, q_AB, ωAB_B) = dq_AB .= q_derivativeCore(q_AB, ωAB_B)


@inline q_rotx!(q, θ) = q .= q_rotx(θ)
@inline q_roty!(q, θ) = q .= q_roty(θ)
@inline q_rotz!(q, θ) = q .= q_rotz(θ)

"""
    q_AB = q_fromAxisAngle(u, θ_AB)

Compute the unitary quaternion given as input an axis-angle representation.
"""
q_fromAxisAngle(u, θ) = iszero(u) ? q_identity() : [cos(0.5θ); sin(0.5θ)*normalize(u)]

function q_fromAxisAngle(u::Vector3, θ)
    iszero(u) && return (1.0, 0.0, 0.0, 0.0)
    x, y, z = normalize(u)
    sθ, cθ = sincos(θ/2)
    return (cθ, x*sθ, y*sθ, z*sθ)
end

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
@inline q_inverse(q::AbstractVector) = collect(q_inverseCore(q))
@inline q_inverse(q::Vector3) = q_inverseCore(q)

@inline q_toRv(q::AbstractVector) = collect(q_toRvCore(q))
@inline q_toRv(q::Vector3)= q_toRvCore(q)

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
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:xyz}) = q_rotx(θ1) * q_roty(θ2) * q_rotz(θ3)#q_fromSequence((θ1, θ2, θ3), (1, 2, 3))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:xzy}) = q_rotx(θ1) * q_rotz(θ2) * q_roty(θ3)#q_fromSequence((θ1, θ2, θ3), (1, 3, 2))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:yxz}) = q_roty(θ1) * q_rotx(θ2) * q_rotz(θ3)#q_fromSequence((θ1, θ2, θ3), (2, 1, 3))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:yzx}) = q_roty(θ1) * q_rotz(θ2) * q_rotx(θ3)#q_fromSequence((θ1, θ2, θ3), (2, 3, 1))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:zxy}) = q_rotz(θ1) * q_rotx(θ2) * q_roty(θ3)#q_fromSequence((θ1, θ2, θ3), (3, 1, 2))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:zyx}) = q_rotz(θ1) * q_roty(θ2) * q_rotx(θ3)#q_fromSequence((θ1, θ2, θ3), (3, 2, 1))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:xyx}) = q_rotx(θ1) * q_roty(θ2) * q_rotx(θ3)#q_fromSequence((θ1, θ2, θ3), (1, 2, 1))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:xzx}) = q_rotx(θ1) * q_rotz(θ2) * q_rotx(θ3)#q_fromSequence((θ1, θ2, θ3), (1, 3, 1))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:yxy}) = q_roty(θ1) * q_rotx(θ2) * q_roty(θ3)#q_fromSequence((θ1, θ2, θ3), (2, 1, 2))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:yzy}) = q_roty(θ1) * q_rotz(θ2) * q_roty(θ3)#q_fromSequence((θ1, θ2, θ3), (2, 3, 2))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:zxz}) = q_rotz(θ1) * q_rotx(θ2) * q_rotz(θ3)#q_fromSequence((θ1, θ2, θ3), (3, 1, 3))
@inline q_fromEuler(θ1, θ2, θ3, ::Val{:zyz}) = q_rotz(θ1) * q_roty(θ2) * q_rotz(θ3)#q_fromSequence((θ1, θ2, θ3), (3, 2, 3))


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

@inline q_attitudeError!(θerr, qNominal, q) = θerr .= q_attitudeErrorCore(qNominal, q)
@inline q_attitudeError(qNominal, q) = q_attitudeErrorCore(qNominal, q)

function q_fromZaxis(zB_A)
    zB_A = normalize(zB_A)
    return q_fromDcm([nullspace(zB_A') zB_A])
end
