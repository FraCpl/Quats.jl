# Author: F. Capolupo
# European Space Agency, 2023
module Quats
using LinearAlgebra

# Useful constants
const ex = (1.0, 0.0, 0.0)
const ey = (0.0, 1.0, 0.0)
const ez = (0.0, 0.0, 1.0)
const I3 = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)

# CAUTION
# Math functions involving Vec3 and Mat3 do not verify the size of the inputs. These shall
# only be used when the user is sure that input data are effectively 3-dimesional vectors
# and/or 3x3 matrices, otherwise the functions may provide an incorrect output.
export Vec3, Quat, Mat3
export toVector, toVector!, toMatrix, toMatrix!
export dot3, cross3, mmul3, mmul3!, mrow3, mcol3, crossSq3, vmul3, vmul3!, getij, evalpoly3
export randn3, zeros3, rand3, randn4, rand4
export crossMat, crossMat!, crossMatSq, crossMatSq!
export cross!, crossSq!, addCross!, addCrossSq!, crossSq
include("math.jl")

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

export dcm_fromAxes,
    dcm_fromAxes!,
    dcm_normalize,
    dcm_rotx,
    dcm_roty,
    dcm_rotz
include("dcm.jl")

export rv_derivative
include("rv.jl")

export transformFrame, transformFrameInv
include("transformFrame.jl")

end
