var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = Quaternions","category":"page"},{"location":"#Quaternions.jl","page":"Home","title":"Quaternions.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package provides a bare-bones representation of attitude transformations through quaternions, rotation vectors, and transformation matrices (a.k.a. direction cosine matrices – DCM). ","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Just type:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(url=\"https://github.com/FraCpl/Quaternions.jl\")","category":"page"},{"location":"#Notation-and-conventions","page":"Home","title":"Notation and conventions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Given a 3-dimensional vector v, a reference frame A and a reference frame B, we indicate with  v^A the projection of v into A, and with v^B its projection in B, so that the two  projections can be related through the transformation matrix R_AB as follows","category":"page"},{"location":"","page":"Home","title":"Home","text":"v^A = R_AB v^B","category":"page"},{"location":"","page":"Home","title":"Home","text":"The same transformation can be represented using the quaternion q_AB","category":"page"},{"location":"","page":"Home","title":"Home","text":"v^A = q_AB  v^B  q_AB^*","category":"page"},{"location":"","page":"Home","title":"Home","text":"where in the formula above, 3D vectors are represented by quaternions having a null scalar component.","category":"page"},{"location":"","page":"Home","title":"Home","text":"In this package, quaternions are represented as 4-dimensional vectors, using the basic vector format natively provided by Julia. The scalar component of the quaternion corresponds to its first element,  i.e., q[1]. We use a right-handed passive convention for representing attitude  transformations through quaternions, and use theta_AB to indicate the angle by which the frame  A needs to be actively rotated to make it coincide with frame B (i.e., the rotation from A to B).","category":"page"},{"location":"","page":"Home","title":"Home","text":"note: Note\nOther (and weirder) quaternions notations exist in the litterature, which are obviously not covered by this package.","category":"page"},{"location":"#API","page":"Home","title":"API","text":"","category":"section"},{"location":"#Quaternions","page":"Home","title":"Quaternions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Modules = [Quaternions]\nOrder   = [:function, :type]\nPages   = [\"q.jl\"]","category":"page"},{"location":"#Quaternions.q_build-Tuple{Any, Vector}","page":"Home","title":"Quaternions.q_build","text":"q_build(qs,qv)\n\nBuild a quaternion from its scalar and vectorial components.\n\n\n\n\n\n","category":"method"},{"location":"#Quaternions.q_derivative-Tuple{Vector, Vector}","page":"Home","title":"Quaternions.q_derivative","text":"q̇_AB = q_derivative(q_AB,ωAB_B)\n\nCompute the time derivative of a unitary quaternion, given the corresponding angular velocity vector.\n\nMathematically, this function performs the following operation:\n\nq_AB = frac12 q_AB  0 ω^B_AB\n\nwhere ωAB_B represents the angular velocity of frame B with respect to frame A, projected into frame B.\n\n\n\n\n\n","category":"method"},{"location":"#Quaternions.q_fromAxes-Tuple{Vector, Vector, Vector}","page":"Home","title":"Quaternions.q_fromAxes","text":"q_AB = q_fromAxes(xB_A, yB_A, zB_A)\n\nCompute the attitude quaternion given as input the axes of a reference frame.\n\n\n\n\n\n","category":"method"},{"location":"#Quaternions.q_fromAxisAngle-Tuple{Vector, Any}","page":"Home","title":"Quaternions.q_fromAxisAngle","text":"q_AB = q_fromAxisAngle(u,θ_AB)\n\nCompute the unitary quaternion given as input an axis-angle representation.\n\n\n\n\n\n","category":"method"},{"location":"#Quaternions.q_fromDcm-Tuple{Matrix}","page":"Home","title":"Quaternions.q_fromDcm","text":"q_AB = q_fromDcm(R_AB)\n\nTranslate the input rotation matrix into a unitary quaternion.\n\n\n\n\n\n","category":"method"},{"location":"#Quaternions.q_identity-Tuple{}","page":"Home","title":"Quaternions.q_identity","text":"q_identity()\n\nGet identity quaternion, with scalar component equal to 1 and vector components equal to zero.\n\n\n\n\n\n","category":"method"},{"location":"#Quaternions.q_inverse-Tuple{Vector}","page":"Home","title":"Quaternions.q_inverse","text":"q⁻¹ = q_inverse(q)\n\nCompute the inverse of the input quaternion.\n\n\n\n\n\n","category":"method"},{"location":"#Quaternions.q_multiply-Tuple{Vector, Vector}","page":"Home","title":"Quaternions.q_multiply","text":"q_AC = q_multiply(q_AB,q_BC)\n\nMutliply the two input quaternions as follows:\n\nq_AC = q_AB  q_BC\n\n\n\n\n\n","category":"method"},{"location":"#Quaternions.q_multiplyn-Tuple","page":"Home","title":"Quaternions.q_multiplyn","text":"q = q_multiplyn(q1,q2,q3,...)\n\n    q = q₁  q₂  q₃  \n\n\n\n\n\n","category":"method"},{"location":"#Quaternions.q_random-Tuple{}","page":"Home","title":"Quaternions.q_random","text":"q_random()\n\nGenerate a random unitary quaternion.\n\n\n\n\n\n","category":"method"},{"location":"#Quaternions.q_toDcm-Tuple{Vector}","page":"Home","title":"Quaternions.q_toDcm","text":"R_AB = q_toDcm(q_AB)\n\nTranslate the input unitary quaternion into a transformation matrix.\n\nR_AB(q_AB) = I + 2qₛqᵥ + 2qᵥ²\n\n\n\n\n\n","category":"method"},{"location":"#Quaternions.q_transformVector-Tuple{Vector, Vector}","page":"Home","title":"Quaternions.q_transformVector","text":"v_A = q_transformVector(q_AB,v_B)\n\n(previously known as q_rotateVector)\n\nProject the vector v from frame B into frame A using the following passive rotation formula\n\nvᴬ = q_AB  vᴮ  q_BA\n\nIn the formula above, 3D vectors are represented by quaternions having a null scalar component.\n\n\n\n\n\n","category":"method"},{"location":"#Quaternions.q_transpose-Tuple{Vector}","page":"Home","title":"Quaternions.q_transpose","text":"q' = q_transpose(q)\n\nTranspose the input quaternion.\n\n\n\n\n\n","category":"method"},{"location":"#Transformation-matrices","page":"Home","title":"Transformation matrices","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Modules = [Quaternions]\nOrder   = [:function, :type]\nPages   = [\"dcm.jl\"]","category":"page"},{"location":"#Quaternions.dcm_fromAxisAngle-Tuple{Vector, Any}","page":"Home","title":"Quaternions.dcm_fromAxisAngle","text":"R_AB = dcm_fromAxisAngle(u,θ_AB)\n\nCompute the transformation matrix given the axis and angle.\n\nR_AB(θ_AB) = I + sin(θ_AB)u + (1 - cos(θ_AB))u^2\n\n\n\n\n\n","category":"method"},{"location":"#Quaternions.dcm_fromQuaternion-Tuple{Vector}","page":"Home","title":"Quaternions.dcm_fromQuaternion","text":"R_AB = dcm_fromQuaternion(q_AB)\n\nCompute a transformation matrix from a quaternion.\n\n\n\n\n\n","category":"method"},{"location":"#Quaternions.dcm_random-Tuple{}","page":"Home","title":"Quaternions.dcm_random","text":"dcm_random()\n\nGenerate a random transformation matrix.\n\n\n\n\n\n","category":"method"},{"location":"#Quaternions.dcm_toQuaternion-Tuple{Matrix}","page":"Home","title":"Quaternions.dcm_toQuaternion","text":"q_AB = dcm_toQuaternion(R_AB)\n\nTranslate a transformation matrix into a quaternion.\n\n\n\n\n\n","category":"method"},{"location":"#Rotation-vectors","page":"Home","title":"Rotation vectors","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Modules = [Quaternions]\nOrder   = [:function, :type]\nPages   = [\"rv.jl\"]","category":"page"},{"location":"#Auxiliary-functions","page":"Home","title":"Auxiliary functions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"crossMat(v::Vector)","category":"page"},{"location":"#Quaternions.crossMat-Tuple{Vector}","page":"Home","title":"Quaternions.crossMat","text":"[v×] = crossMat(v)\n\nCompute the cross product matrix of a vector v.\n\n\n\n\n\n","category":"method"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"}]
}
