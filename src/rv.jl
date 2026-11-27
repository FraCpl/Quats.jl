@inline function rv_derivative(ϕ_AB, ωAB_B)
    # Bortz equation
    θ = norm(ϕ_AB)
    ϕxω = cross3(ϕ_AB, ωAB_B)
    k = 1/θ^2 - 1 / (2 * θ * tan(θ))
    return ωAB_B + ϕxω ./ 2 + k .* cross3(ϕ_AB, ϕxω)
end
