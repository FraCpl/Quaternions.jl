"""
    dcm_random()

Generate a random transformation matrix.
"""
dcm_random() = dcm_fromAxisAngle(randn(3), 2π*rand())

"""
    R_AB = dcm_fromAxisAngle(u, θ_AB)

Compute the transformation matrix given the axis and angle.
```math
R_{AB}(θ_{AB}) = I + \\sin(θ_{AB})[u×] + (1 - \\cos(θ_{AB}))[u×]^2
```
"""
function dcm_fromAxisAngle(u, θ)
    ux = crossMat(normalize(u))
    return I + sin(θ)*ux + (1 - cos(θ))*ux*ux
end

dcm_fromAxisAngle(idx::Int, θ) = dcm_fromAxisAngle(Float64.([idx==1; idx==2; idx==3]), θ)

# θ = [θ_AB, θ_BC, θ_CD] --> R_AD
function dcm_fromEuler(sequence::Vector{Int}, θ)
    R = I
    for i in eachindex(sequence)
        R = R*dcm_fromAxisAngle(sequence[i], θ[i])
    end
    return R
end

"""
    q_AB = dcm_toQuaternion(R_AB)

Translate a transformation matrix into a quaternion.
"""
dcm_toQuaternion(R::Matrix) = q_fromDcm(R)

"""
    R_AB = dcm_fromQuaternion(q_AB)

Compute a transformation matrix from a quaternion.
"""
dcm_fromQuaternion(q) = q_toDcm(q)

dcm_fromRv(ϕ) = dcm_fromAxisAngle(ϕ, norm(ϕ))
dcm_toRv(R::Matrix) = q_toRv(dcm_toQuaternion(R))
