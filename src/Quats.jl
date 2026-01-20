# Author: F. Capolupo
# European Space Agency, 2023
module Quats
using LinearAlgebra
# using StaticArrays

export crossMat, crossMat!, crossMatInv, crossMatSq, crossMatSq!
export cross!, crossSq!, addCross!, addCrossSq!, crossMatStatic, crossMatSqStatic
include("cross.jl")

export q_multiply,
    q_fromAxes,
    q_random,
    q_transformVector,
    q_transpose,
    q_fromDcm,
    q_derivative,
    q_multiplyn,
    q_fromAxisAngle,
    q_fromAxisAngle!,
    q_toAxes,
    q_toDcm,
    q_inverse,
    q_build,
    q_identity,
    q_toRv,
    q_fromRv,
    q_fromRv!,
    q_rate,
    q_toEuler,
    q_fromEuler,
    q_interp,
    q_slerp,
    q_toAxisAngle,
    q_attitudeError,
    q_testConvention,
    q_derivative!,
    q_transpose!,
    q_transformVector!,
    q_multiply!,
    q_toDcm!,
    q_fromDcm!,
    q_attitudeError!,
    q_tox,
    q_tox!,
    q_toy,
    q_toy!,
    q_toz,
    q_toz!,
    q_fromAxes!,
    q_transformVectorT!,
    q_multiplyT1!,
    q_multiplyT2!,
    q_multiplyT12!
include("q.jl")
#include("q_static.jl")

export dcm_random,
    dcm_fromAxisAngle,
    dcm_toQuaternion,
    dcm_fromQuaternion,
    dcm_toEuler,
    dcm_toRv,
    dcm_fromRv,
    dcm_fromEuler,
    dcm_fromAxes,
    dcm_rotAxis,
    dcm_normalize
include("dcm.jl")

export rv_toQuaternion, rv_fromQuaternion, rv_derivative, rv_fromDcm, rv_toDcm
include("rv.jl")

end
