rv_toQuaternion(ϕ::Vector) = q_fromRv(ϕ)
rv_fromQuaternion(q::Vector) = q_toRv(q)
rv_toDcm(ϕ::Vector) = dcm_fromRv(ϕ)
rv_fromDcm(R::Matrix) = dcm_toRv(R)

function rv_derivative(ϕ_AB::Vector, ωAB_B::Vector)
    # Bortz equation
    θ = norm(ϕ_AB)
    ϕxω = ϕ × ωAB_B
    return ωAB_B + 0.5*ϕxω + (1.0/θ^2 - 0.5/θ/tan(θ))* (ϕ × ϕxω)
end
