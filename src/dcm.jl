"""
    dcm_random()

Generate a random transformation matrix.
"""
dcm_random() = q_toDcm(q_random())

"""
    R_AB = dcm_fromAxisAngle(u, θ_AB)

Compute the transformation matrix given the axis and angle.
```math
R_{AB}(θ_{AB}) = I + \\sin(θ_{AB})[u×] + (1 - \\cos(θ_{AB}))[u×]^2
```
"""
dcm_fromAxisAngle(u, θ) = q_toDcm(q_fromAxisAngle(u, θ))
dcm_fromAxisAngle(idx::Int, θ) = q_toDcm(q_fromAxisAngle(idx, θ))

# θ = [θ_AB, θ_BC, θ_CD] --> R_AD
dcm_fromEuler(θ, sequence::Vector{Int}=[3, 2, 1]) = q_toDcm(q_fromEuler(θ, sequence))
dcm_toEuler(R, sequence::Vector{Int}=[3, 2, 1]) = q_toEuler(q_fromDcm(R), sequence)

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

function dcm_fromRv(ϕ)
    θ = norm(ϕ)
    if θ == 0.0
        return Matrix(1.0I, 3, 3)
    end
    return dcm_fromAxisAngle(ϕ, θ)
end

dcm_toRv(R::Matrix) = q_toRv(dcm_toQuaternion(R))

"""
    R_AB = dcm_fromAxes(xB_A, yB_A, zB_A)

Compute the transformation matrix given as input the axes of a reference frame.
"""
function dcm_fromAxes(xB_A, yB_A, zB_A)
    if isempty(xB_A); xB_A = yB_A × zB_A; end
    if isempty(yB_A); yB_A = zB_A × xB_A; end
    if isempty(zB_A); zB_A = xB_A × yB_A; end

    return [xB_A yB_A zB_A] # R_AB
end
