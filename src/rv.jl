rv_toQuaternion(ϕ) = q_fromRv(ϕ)
rv_fromQuaternion(q) = q_toRv(q)
rv_toDcm(ϕ) = dcm_fromRv(ϕ)
rv_fromDcm(R) = dcm_toRv(R)

function rv_derivative(ϕ_AB, ωAB_B)
    # Bortz equation
    θ = norm(ϕ_AB)
    ϕxω = ϕ_AB × ωAB_B
    return ωAB_B + 0.5*ϕxω + (1.0/θ^2 - 0.5/θ/tan(θ))* (ϕ_AB × ϕxω)
end
