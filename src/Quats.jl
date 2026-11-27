# Author: F. Capolupo
# European Space Agency, 2023
module Quats
using LinearAlgebra

# Useful constants
const ex = (1.0, 0.0, 0.0)
const ey = (0.0, 1.0, 0.0)
const ez = (0.0, 0.0, 1.0)
const I3 = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)

# Core types definition
const Vec3{T} = NTuple{3, T}        # 3x1 Vector
const Quat{T} = NTuple{4, T}        # 4x1 Quaternion
const Mat3{T} = NTuple{9, T}        # 3x3 Matrix
export Vec3, Quat, Mat3

# CAUTION 1:
# This package is "contagious", i.e., when hybrid inputs are provided, the output is always
# one of the custom types defined by the package, i.e., Vec3, Mat3, or Quat.

# CAUTION 2:
# Math functions involving Vec3 and Mat3 do not verify the size of the inputs. These shall
# only be used when the user is sure that input data are effectively 3-dimesional vectors
# and/or 3x3 matrices, otherwise the functions may provide an incorrect output.
export toVector, toMatrix,
    mrow3, mcol3, vmul3, vmul3!, getij, evalpoly3,
    zeros3, randn3, rand3, randn4, rand4,
    cross!, addCross!, crossSq, crossSq!, addCrossSq!,
    crossMat, crossMat!, crossMatSq, crossMatSq!
include("vec3.jl")
include("mat3.jl")

export q_zero, q_random, toVector,
    q_multiply, q_multiply!, q_multiplyT1!, q_multiplyT2!, q_multiplyT12!, q_multiplyn, q_multiplyn!,
    q_toDcm, q_toDcm!,
    q_fromDcm, q_fromDcm!,
    q_transformVector, q_transformVector!, q_transformVectorT, q_transformVectorT!
    q_transpose, q_transpose!,
    q_derivative, q_derivative!,
    q_rotx, q_rotx!, q_roty, q_roty!, q_rotz, q_rotz!,
    q_fromAxisAngle, q_fromAxisAngle!,
    q_toAxisAngle, q_fromAxes, q_fromAxes!, q_toAxes,
    q_inverse,
    q_toRv,
    q_fromRv, q_fromRv!,
    q_rate,
    q_toEuler, q_fromEuler,
    q_interp, q_slerp,
    q_attitudeError, q_attitudeError!,
    q_tox, q_toy, q_toz
include("quat.jl")

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
