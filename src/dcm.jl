"""
    R_AB = dcm_fromAxes(xB_A, yB_A, zB_A)

Compute the transformation matrix given as input the axes of a reference frame.
"""
@inline dcm_fromAxes(xB_A::AbstractVector, yB_A::AbstractVector, zB_A::AbstractVector) = [xB_A yB_A zB_A] # R_AB
@inline dcm_fromAxes(xB_A, yB_A, zB_A) = (xB_A[1], xB_A[2], xB_A[3], yB_A[1], yB_A[2], yB_A[3], zB_A[1], zB_A[2], zB_A[3])

@inline function dcm_fromAxes!(R_AB, xB_A, yB_A, zB_A)
    @inbounds for i in 1:3
        R_AB[i, 1] = xB_A[i]
        R_AB[i, 2] = yB_A[i]
        R_AB[i, 3] = zB_A[i]
    end
    return R_AB
end

function dcm_rotx(angle)
    sθ, cθ = sincos(angle)
    return (1.0, 0.0, 0.0, 0.0, cθ, sθ, 0.0, -sθ, cθ)
end

function dcm_roty(angle)
    sθ, cθ = sincos(angle)
    return (cθ, 0.0, -sθ, 0.0, 1.0, 0.0, sθ, 0.0, cθ)
end

function dcm_rotz(angle)
    sθ, cθ = sincos(angle)
    return (cθ, sθ, 0.0, -sθ, cθ, 0.0, 0.0, 0.0, 1.0)
end

@inline function dcm_normalize(R::AbstractMatrix)
    U, ~, V = svd(R; full=true)
    return U * V'
end

@inline dcm_normalize(R::Mat3) = Mat3(dcm_normalize(toMatrix(R)))
