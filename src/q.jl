"""
    q_AC = q_multiply(q_AB, q_BC)

Mutliply the two input quaternions as follows:
```math
q_{AC} = q_{AB} ⊗ q_{BC}
```
"""
function q_multiply(q_AB, q_BC)
    ps = q_AB[1]; pv = q_AB[2:4]
    qs = q_BC[1]; qv = q_BC[2:4]
    return [ps*qs - (pv ⋅ qv); ps.*qv + qs.*pv + pv × qv]  # = q_AC = p x q, p = q_AB, q = q_BC
end

"""
    q = q_multiplyn(q1, q2, q3, ...)

```math
    q = q₁ ⊗ q₂ ⊗ q₃ ⊗ ...
```
"""
function q_multiplyn(q...)
    qOut = q[1]
    for i in 2:lastindex(q)
        qOut = q_multiply(qOut, q[i])
    end

    return qOut
end

"""
    q_build(qs,qv)

Build a quaternion from its scalar and vectorial components.
"""
q_build(qs, qv) = [qs; qv]

"""
    R_AB = q_toDcm(q_AB)

Translate the input unitary quaternion into a transformation matrix.
```math
R_{AB}(q_{AB}) = I + 2qₛ[qᵥ×] + 2[qᵥ×]²
```
"""
function q_toDcm(q)     # R_BA from q_BA
    qx = crossMat(q[2:4])
    return I + 2.0*(qx*qx + q[1].*qx)
end

"""
    q_AB = q_fromDcm(R_AB)

Translate the input rotation matrix into a unitary quaternion.
"""
function q_fromDcm(R_BA)
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
q_fromAxes(xB_A, yB_A, zB_A) = q_fromDcm(dcm_fromAxes(xB_A, yB_A, zB_A))

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
function q_transformVector(q_AB, v_B)
    qxv = q_AB[2:4] × v_B
    return v_B + 2.0*(q_AB[2:4] × qxv + q_AB[1].*qxv) # v_A
end

"""
    q' = q_transpose(q)

Transpose the input quaternion.
"""
q_transpose(q) = [q[1]; -q[2:4]]

"""
    q̇_AB = q_derivative(q_AB, ωAB_B)

Compute the time derivative of a unitary quaternion, given the corresponding
angular velocity vector.

Mathematically, this function performs the following operation:
```math
q̇_{AB} = \\frac{1}{2} q_{AB} ⊗ [0; ω^B_{AB}]
```
where ```ωAB_B``` represents the angular velocity of frame ``B`` with respect to
frame ``A``, projected into frame ``B``.
"""
q_derivative(q_AB, ωAB_B) = q_multiply(q_AB, [0.0; 0.5.*ωAB_B])  # dq_BA

"""
    q_AB = q_fromAxisAngle(u, θ_AB)

Compute the unitary quaternion given as input an axis-angle representation.
"""
q_fromAxisAngle(u, θ) = [cos(0.5θ); sin(0.5θ)*normalize(u)]
q_fromAxisAngle(idx::Int, θ) = q_fromAxisAngle(Float64.([idx==1; idx==2; idx==3]), θ)

function q_toAxisAngle(q)
    nqv = norm(q[2:4])
    θ = 2atan(nqv, q[1])
    return q[2:4]./nqv, θ
end

function q_toAxes(q_AB)
    q_BA = q_transpose(q_AB)
    xB_A = q_transformVector(q_BA,[1.0; 0.0; 0.0])
    yB_A = q_transformVector(q_BA,[0.0; 1.0; 0.0])
    zB_A = q_transformVector(q_BA,[0.0; 0.0; 1.0])

    return xB_A, yB_A, zB_A
end

"""
    q⁻¹ = q_inverse(q)

Compute the inverse of the input quaternion.
"""
q_inverse(q) = q_transpose(q)./(q ⋅ q)


function q_toRv(q)
    nqv = norm(q[2:4])
    if nqv < 1e-10
        return zeros(3)
    end
    return (2*atan(nqv, q[1])/nqv).*q[2:4]
end

q_fromRv(ϕ) = q_fromAxisAngle(ϕ, norm(ϕ))

"""
This uses rotation vector to compute the average angular rate between two
sampled quaternion values, knowing that:
q[k] = q[k-1] o dq
and
angRate(t) = dtheta/dt for t in [t[k-1],t[k]]
This means that angRate is the equivalent constant average rate that
rotates q[k-1] into q[k].
"""
@views function q_rate(t, q_AB)
    dt = diff(t)
    dq_AB = q_multiply.(q_transpose.(q_AB[1:end-1]), q_AB[2:end])
    return [q_toRv.(dq_AB)./dt; [zeros(3)]]    # angRateAB_B
end

# https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
# 3(yaw)-2(pitch)-1(roll) sequence
# CAUTION: Only [3, 2, 1] and [1, 2, 3] sequencs are implemented for now
function q_toEuler(q, sequence=[3, 2, 1])
    qs, qx, qy, qz = q

    if sequence == [3, 2, 1]
        # roll (x-axis rotation)
        t1 = 2(qs*qx + qy*qz)
        t2 = 1 - 2(qx*qx + qy*qy)
        roll = atan(t1, t2)

        # pitch (y-axis rotation)
        t3 = 2(qs*qy - qz*qx)
        t3 = max(-1.0, min(t3, 1.0))
        pitch = asin(t3)

        # yaw (z-axis rotation)
        t4 = 2(qs*qz + qx*qy)
        t5 = 1 - 2(qy*qy + qz*qz)
        yaw = atan(t4, t5)

        return [yaw; pitch; roll]

    elseif sequence == [1, 2, 3]
        sinPitch = 2(qx*qz + qy*qs)
        mask = sinPitch ≥ 1 - 10*eps() || sinPitch ≤ -1 + 10*eps()
        sinPitch = max(-1.0, min(sinPitch, 1.0))
        pitch = asin(sinPitch)
        if mask
            roll = 0
            yaw = sign(sinPitch)*2*atan(qx, qs)
        else
            roll = atan(-2(qy*qz - qx*qs), qs^2 - qx^2 - qy^2 + qz^2);
            yaw = atan(-2(qx*qy - qz*qs), qs^2 + qx^2 - qy^2 - qz^2)
        end
        return [roll; pitch; yaw]
    end
end

function q_fromEuler(θ, sequence=[3, 2, 1])
    q = q_fromAxisAngle(sequence[1], θ[1])
    for i in 2:lastindex(θ)
        q .= q_multiply(q, q_fromAxisAngle(sequence[i], θ[i]))
    end
    return q
end

function q_exp(q)
    qvn = norm(q[2:4])
    qvnn = qvn
    if qvnn == 0.0
        qvnn = 1.0
    end
    qv = q[2:4]./qvnn
    return exp(q[1]).*[cos(qvn); qv.*sin(qvn)]
end

function q_log(q)
    vNorm = norm(q[2:4])
    if vNorm == 0.0
        vNorm = 1.0
    end
    nq = norm(q)
    return [log(nq); q[2:4]/vNorm.*acos(q[1]/nq)]
end

q_power(q, n) = q_exp(n.*q_log(q))

# τ in [0, 1]
q_slerp(q0, q1, τ) = q_multiply(q0, q_power(q_multiply(q_transpose(q0), q1), τ))

function q_interp(t, q, ti)
    id0 = findlast(ti .≥ t)
    id1 = min(id0+1, length(t))
    return q_slerp(q[id0], q[id1], (ti - t[id0])/(t[id1] - t[id0]))
end

# e.g., qNominal = qDes, q = qEst
function q_attitudeError(qNominal, q)
    imax = findmax(abs.(qNominal))[2]
    sgn = sign(qNominal[imax]) == sign(q[imax]) ? 1.0 : -1.0
    return 2q_multiply(q_transpose(sgn*q), qNominal)[2:4]
end


# This function determines the convention used by a quaternion multiplication function in
# terms of: 1) Scalar-first or scalar-last, 2) Right-handed or left-handed algebra
# The multiplication function shall provide q*p = fmult(q, p)
function q_testConvention(fmult)
    q = [1.0; 2.0; 3.0; 4.0]
    p = [-4.0; 2.0; -3.0; 1.0]
    q_testConvention(fmult(q, p), q, p)
end

function q_testConvention(qp, q, p)

    # out = q*p
    function q_multiplyUniversal(q, p; scalarFirst=true, right=true)
        qq = copy(q); pp = copy(p)
        if !scalarFirst
            qq = [q[4]; q[1:3]]
            pp = [p[4]; p[1:3]]
        end
        qqpp = right ? qp_right(qq, pp) : qp_left(qq, pp)
        if scalarFirst; return qqpp; end
        return [qqpp[2:4]; qqpp[1]]
    end

    function qp_right(q, p)
        qs, qx, qy, qz = q
        return [+qs -qx -qy -qz;
                +qx +qs -qz +qy;
                +qy +qz +qs -qx;
                +qz -qy +qx +qs]*p
    end

    function qp_left(q, p)
        qs, qx, qy, qz = q
        return [+qs -qx -qy -qz;
                +qx +qs +qz -qy;
                +qy -qz +qs +qx;
                +qz +qy -qx +qs]*p
    end

    qerr(a, b) = min(norm(a - b), norm(a + b))
    qp1 = q_multiplyUniversal(q, p; scalarFirst=true, right=true)
    qp2 = q_multiplyUniversal(q, p; scalarFirst=false, right=true)
    qp3 = q_multiplyUniversal(q, p; scalarFirst=true, right=false)
    qp4 = q_multiplyUniversal(q, p; scalarFirst=false, right=false)
    errv = [qerr(qp1, qp), qerr(qp2, qp), qerr(qp3, qp), qerr(qp4, qp)]

    msgs = ["Scalar first, right-handed (ij = k)", "Scalar last, right-handed (ij = k)", "Scalar first, left-handed (ij = -k)", "Scalar-last, left-handed (ij = -k)"]
    println("Quaternion convention: \n"*msgs[findmin(errv)[2]])
end
