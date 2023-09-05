"""
    q_AC = q_multiply(q_AB,q_BC)

Mutliply the two input quaternions as follows:
```math
q_{AC} = q_{AB} ⊗ q_{BC}
```
"""
function q_multiply(q_AB::Vector, q_BC::Vector)
    ps = q_AB[1]; pv = q_AB[2:4]
    qs = q_BC[1]; qv = q_BC[2:4]
    return [ps*qs - (pv ⋅ qv); ps.*qv + qs.*pv + pv × qv]  # = q_AC = p x q, p = q_AB, q = q_BC
end

"""
    q = q_multiplyn(q1,q2,q3,...)

```math
    q = q₁ ⊗ q₂ ⊗ q₃ ⊗ ...
```
"""
function q_multiplyn(q...)
    qOut = q[1]
    for i in 2:lastindex(q)
        qOut = q_multiply(qOut,q[i])
    end

    return qOut
end

"""
    q_build(qs,qv)

Build a quaternion from its scalar and vectorial components.
"""
q_build(qs,qv::Vector) = [qs; qv]

"""
    R_AB = q_toDcm(q_AB)

Translate the input unitary quaternion into a transformation matrix.
```math
R_{AB}(q_{AB}) = I + 2qₛ[qᵥ×] + 2[qᵥ×]²
```
"""
function q_toDcm(q::Vector)     # R_BA from q_BA
    qx = crossMat(q[2:4])
    return I + 2.0*(qx*qx + q[1].*qx)
end

"""
    q_AB = q_fromDcm(R_AB)

Translate the input rotation matrix into a unitary quaternion.
"""
function q_fromDcm(R_BA::Matrix)
    dcm11 = R_BA[1,1]; dcm12 = R_BA[2,1]; dcm13 = R_BA[3,1];
    dcm21 = R_BA[1,2]; dcm22 = R_BA[2,2]; dcm23 = R_BA[3,2];
    dcm31 = R_BA[1,3]; dcm32 = R_BA[2,3]; dcm33 = R_BA[3,3];

    vv = 1.0 .+ [+dcm11-dcm22-dcm33; -dcm11+dcm22-dcm33; -dcm11-dcm22+dcm33; +dcm11+dcm22+dcm33]
    idx = argmax(vv);
    qx = 0.5*sqrt(abs(vv[idx]))
    f  = 0.25/qx;

    if idx == 1
        return [f*(dcm23 - dcm32); qx; f*(dcm12 + dcm21);  f*(dcm31 + dcm13)]
    elseif idx == 2
        return [f*(dcm31 - dcm13); f*(dcm12 + dcm21); qx; f*(dcm23 + dcm32)]
    elseif idx == 3
        return [f*(dcm12 - dcm21); f*(dcm31 + dcm13); f*(dcm23 + dcm32); qx]
    end
    return [qx; f*(dcm23 - dcm32); f*(dcm31 - dcm13); f*(dcm12 - dcm21)]
end

"""
    q_AB = q_fromAxes(xB_A, yB_A, zB_A)

Compute the attitude quaternion given as input the axes of a reference frame.
"""
function q_fromAxes(xB_A::Vector, yB_A::Vector, zB_A::Vector)
    if isempty(xB_A); xB_A = yB_A × zB_A; end
    if isempty(yB_A); yB_A = zB_A × xB_A; end
    if isempty(zB_A); zB_A = xB_A × yB_A; end

    return q_fromDcm([xB_A yB_A zB_A])
end

"""
    q_random()

Generate a random unitary quaternion.
"""
q_random() = normalize(randn(4))

"""
    q_identity()

Get identity quaternion, with scalar component equal to 1 and
vector components equal to zero.
"""
q_identity() = [1.0; 0.0; 0.0; 0.0]


"""
    v_A = q_transformVector(q_AB,v_B)
(previously known as q_rotateVector)

Project the vector v from frame B into frame A using the following
passive rotation formula
```math
vᴬ = q_{AB} ⊗ vᴮ ⊗ q_{BA}
```
In the formula above, 3D vectors are represented by quaternions having
a null scalar component.
"""
function q_transformVector(q_AB::Vector, v_B::Vector)
    qxv = q_AB[2:4] × v_B
    return v_B + 2.0*(q_AB[2:4] × qxv + q_AB[1].*qxv) # v_A
end

"""
    q' = q_transpose(q)

Transpose the input quaternion.
"""
q_transpose(q::Vector) = [q[1]; -q[2:4]]

"""
    q̇_AB = q_derivative(q_AB,ωAB_B)

Compute the time derivative of a unitary quaternion, given the corresponding
angular velocity vector.

Mathematically, this function performs the following operation:
```math
q̇_{AB} = \\frac{1}{2} q_{AB} ⊗ [0; ω^B_{AB}]
```
where ```ωAB_B``` represents the angular velocity of frame ``B`` with respect to
frame ``A``, projected into frame ``B``.
"""
q_derivative(q_AB::Vector, ωAB_B::Vector) = q_multiply(q_AB,[0.0; 0.5.*ωAB_B])  # dq_BA

"""
    q_AB = q_fromAxisAngle(u,θ_AB)

Compute the unitary quaternion given as input an axis-angle representation.
"""
q_fromAxisAngle(u::Vector, θ) = [cos(0.5θ); sin(0.5θ)*normalize(u)]

q_fromAxisAngle(idx::Int, θ) = q_fromAxisAngle(Float64.([idx==1; idx==2; idx==3]), θ)

function q_toAxes(q_BA::Vector)
    xB_A = q_rotateVector(q_BA,[1.0; 0.0; 0.0])
    yB_A = q_rotateVector(q_BA,[0.0; 1.0; 0.0])
    zB_A = q_rotateVector(q_BA,[0.0; 0.0; 1.0])

    return xB_A, yB_A, zB_A
end

"""
    q⁻¹ = q_inverse(q)

Compute the inverse of the input quaternion.
"""
q_inverse(q::Vector) = q_transpose(q)./(q ⋅ q)