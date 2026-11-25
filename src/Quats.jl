# Author: F. Capolupo
# European Space Agency, 2023
module Quats
using LinearAlgebra

const ex = (1.0, 0.0, 0.0)
const ey = (0.0, 1.0, 0.0)
const ez = (0.0, 0.0, 1.0)
const I3 = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)

export Vector3, Quaternion, Matrix3x3
export toVector, toVector!, toMatrix, toMatrix!
export dot3, cross3, mmul3, mmul3!, mrow3, mcol3, crossSq3, dprod3, getij
export randn3, zeros3, rand3, evalpoly3
include("math.jl")

export crossMat, crossMat!, crossMatSq, crossMatSq!
export cross!, crossSq!, addCross!, addCrossSq!, crossMatStatic, crossMatSqStatic, crossSq
include("cross.jl")

include("qCore.jl")
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
    q_multiplyn!,
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
    q_transformVectorT,
    q_multiplyT1!,
    q_multiplyT2!,
    q_multiplyT12!,
    q_fromSequence,
    q_rotx,
    q_roty,
    q_rotz,
    q_rotx!,
    q_roty!,
    q_rotz!
include("q.jl")

export dcm_random,
    dcm_fromAxisAngle,
    dcm_toQuaternion,
    dcm_fromQuaternion,
    dcm_toEuler,
    dcm_toRv,
    dcm_fromRv,
    dcm_fromEuler,
    dcm_fromAxes,
    dcm_fromAxes!,
    dcm_rotAxis,
    dcm_normalize,
    dcm_rotx,
    dcm_roty,
    dcm_rotz
include("dcm.jl")

export rv_toQuaternion, rv_fromQuaternion, rv_derivative, rv_fromDcm, rv_toDcm
include("rv.jl")

export transformFrame, transformFrameInv
include("transformFrame.jl")

end
