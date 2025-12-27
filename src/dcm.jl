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
@inline dcm_fromAxisAngle(idx::Int, θ) = dcm_fromEuler([θ], [idx])

# θ = [θ_AB, θ_BC, θ_CD] --> R_AD
function dcm_fromEuler(θ, sequence::Vector{Int}=[3, 2, 1])
    R = dcm_rotAxis(θ[1], sequence[1])
    for k in 2:lastindex(sequence)
        R .= R*dcm_rotAxis(θ[k], sequence[k])
    end
    return R
end
@inline dcm_toEuler(R, sequence::Vector{Int}=[3, 2, 1]) = q_toEuler(q_fromDcm(R), sequence)

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
@inline function dcm_fromAxes(xB_A, yB_A, zB_A)
    if isempty(xB_A)
        ;
        xB_A = yB_A × zB_A;
    end
    if isempty(yB_A)
        ;
        yB_A = zB_A × xB_A;
    end
    if isempty(zB_A)
        ;
        zB_A = xB_A × yB_A;
    end

    return [xB_A yB_A zB_A] # R_AB
end

# R_AB(θ_AB)
@inline function dcm_rotAxis(angle, axis::Int)
    sθ, cθ = sincos(angle)
    if axis == 1
        return [1 0 0; 0 cθ -sθ; 0 sθ cθ]
    elseif axis == 2
        return [cθ 0 sθ; 0 1 0; -sθ 0 cθ]
    end
    return [cθ -sθ 0; sθ cθ 0; 0 0 1]
end

@inline function dcm_normalize(R)
    U, ~, V = svd(R; full=true)
    return U*V'
end
