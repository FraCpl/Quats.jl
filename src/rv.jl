@inline rv_toQuaternion(ϕ) = q_fromRv(ϕ)
@inline rv_fromQuaternion(q) = q_toRv(q)
@inline rv_toDcm(ϕ) = dcm_fromRv(ϕ)
@inline rv_fromDcm(R) = dcm_toRv(R)

@inline function rv_derivative(ϕ_AB, ωAB_B)
    # Bortz equation
    θ = norm(ϕ_AB)
    ϕxω = ϕ_AB × ωAB_B
    k = 1/θ^2 - 0.5/θ/tan(θ)
    return ωAB_B + 0.5 .* ϕxω + k .* (ϕ_AB × ϕxω)
end
