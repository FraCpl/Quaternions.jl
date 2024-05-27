# Author: F. Capolupo
# European Space Agency, 2023
module Quaternions
using LinearAlgebra

export crossMat
"""
    [v×] = crossMat(v)

Compute the cross product matrix of a vector ```v```.
"""
crossMat(v) = [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]

export crossMatInv
crossMatInv(M::Matrix) = [-M[2, 3]; M[1, 3]; -M[1, 2]]

export q_multiply, q_fromAxes, q_random, q_transformVector,
    q_transpose, q_fromDcm, q_derivative, q_multiplyn,
    q_fromAxisAngle, q_toAxes, q_toDcm, q_inverse, q_build,
    q_identity, q_toRv, q_fromRv, q_rate, q_toEuler, q_fromEuler,
    q_interp, q_slerp, q_toAxisAngle, q_attitudeError
include("q.jl")

export dcm_random, dcm_fromAxisAngle, dcm_toQuaternion, dcm_fromQuaternion, dcm_toEuler,
    dcm_toRv, dcm_fromRv, dcm_fromEuler, dcm_fromAxes
include("dcm.jl")

export rv_toQuaternion, rv_fromQuaternion, rv_derivative, rv_fromDcm, rv_toDcm
include("rv.jl")

end
