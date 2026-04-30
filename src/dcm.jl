"""
    dcm_random()

Generate a random transformation matrix.
"""
@inline dcm_random() = q_toDcm(q_random())

"""
    R_AB = dcm_fromAxisAngle(u, θ_AB)

Compute the transformation matrix given the axis and angle.
```math
R_{AB}(θ_{AB}) = I + \\sin(θ_{AB})[u×] + (1 - \\cos(θ_{AB}))[u×]^2
```
"""
@inline dcm_fromAxisAngle(u, θ) = q_toDcm(q_fromAxisAngle(u, θ))
@inline dcm_fromAxisAngle(idx::Int, θ) = dcm_fromSequence((θ,), (idx,))

# (θ1=θ_AB, θ2=θ_BC, θ3=θ_CD) --> R_AD
@inline dcm_fromEuler(θ1, θ2, θ3, s::Symbol=:zyx) = q_toDcm(q_fromEuler(θ1, θ2, θ3, s))
@inline dcm_fromSequence(θ, sequence=(3, 2, 1)) = q_toDcm(q_fromSequence(θ, sequence))
@inline dcm_toEuler(R, s::Symbol=:zyx) = q_toEuler(q_fromDcm(R), s)

"""
    q_AB = dcm_toQuaternion(R_AB)

Translate a transformation matrix into a quaternion.
"""
@inline dcm_toQuaternion(R::Matrix) = q_fromDcm(R)

"""
    R_AB = dcm_fromQuaternion(q_AB)

Compute a transformation matrix from a quaternion.
"""
@inline dcm_fromQuaternion(q) = q_toDcm(q)

function dcm_fromRv(ϕ)
    θ = norm(ϕ)
    if θ == 0.0
        return Matrix(1.0I, 3, 3)
    end
    return dcm_fromAxisAngle(ϕ, θ)
end

@inline dcm_toRv(R::Matrix) = q_toRv(dcm_toQuaternion(R))

"""
    R_AB = dcm_fromAxes(xB_A, yB_A, zB_A)

Compute the transformation matrix given as input the axes of a reference frame.
"""
@inline dcm_fromAxes(xB_A, yB_A, zB_A) = [xB_A yB_A zB_A] # R_AB

@inline function dcm_fromAxes!(R_AB, xB_A, yB_A, zB_A)
    @inbounds for i in 1:3
        R_AB[i, 1] = xB_A[i]
        R_AB[i, 2] = yB_A[i]
        R_AB[i, 3] = zB_A[i]
    end
    return R_AB
end

# R_AB(θ_AB)
@inline function dcm_rotAxis(angle, axis::Int)
    if axis == 1
        return dcm_rotx(angle)
    elseif axis == 2
        return dcm_roty(angle)
    end
    return dcm_rotz(angle)
end

function dcm_rotx(angle)
    sθ, cθ = sincos(angle)
    return [1 0 0; 0 cθ -sθ; 0 sθ cθ]
end

function dcm_roty(angle)
    sθ, cθ = sincos(angle)
    return [cθ 0 sθ; 0 1 0; -sθ 0 cθ]
end

function dcm_rotz(angle)
    sθ, cθ = sincos(angle)
    return [cθ -sθ 0; sθ cθ 0; 0 0 1]
end

@inline function dcm_normalize(R)
    U, ~, V = svd(R; full=true)
    return U*V'
end
