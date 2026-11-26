# TYPES DEFINITIONS ---------------------------------------------------------------------- #
const Vector3{T} = NTuple{3, T}
const Quaternion{T} = NTuple{4, T}
const Matrix3x3{T} = NTuple{9, T}

# CORE MATH FUNCTIONS -------------------------------------------------------------------- #
@inline zeros3() = (0.0, 0.0, 0.0)
@inline randn3() = (randn(), randn(), randn())
@inline randn3(rng) = (randn(rng), randn(rng), randn(rng))
@inline rand3() = (rand(), rand(), rand())
@inline rand3(rng) = (rand(rng), rand(rng), rand(rng))
@inline randn4() = (randn(), randn(), randn(), randn())
@inline randn4(rng) = (randn(rng), randn(rng), randn(rng), randn(rng))
@inline sum3(a, b) = a[1]+b[1], a[2]+b[2], a[3]+b[3]
@inline subt3(a, b) = a[1]-b[1], a[2]-b[2], a[3]-b[3]
@inline smult3(s, a) = s*a[1], s*a[2], s*a[3]
@inline dot3(a, b) = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
@inline cross3(a, b) = a[2]*b[3] - a[3]*b[2], a[3]*b[1] - a[1]*b[3], a[1]*b[2] - a[2]*b[1]       # out = a × (a × b)

# R * v where R is a 3x3 matrix, and v is a 3x1 vector
@inline function mmul3(R, v)
    r11, r21, r31, r12, r22, r32, r13, r23, r33 = R
    x, y, z = v
    ox = r11 * x + r12 * y + r13 * z
    oy = r21 * x + r22 * y + r23 * z
    oz = r31 * x + r32 * y + r33 * z
    return ox, oy, oz
end

@inline mmul3!(y, R, v) = y .= mmul3(R, v)
@inline dprod3core(a, b) = (a[1]*b[1], a[2]*b[1], a[3]*b[1], a[1]*b[2], a[2]*b[2], a[3]*b[2], a[1]*b[3], a[2]*b[3], a[3]*b[3])

# Dyadic product of two vectors, i.e., a * b'
@inline dprod3(a::AbstractVector, b::AbstractVector=a) = toMatrix(dprod3core(a, b))
@inline dprod3(a, b=a) = dprod3core(a, b)

function crossSq3(a, b)
    a1, a2, a3 = a
    b1, b2, b3 = b
    o1 = (-a2*a2 - a3*a3)*b1 + a1*a2*b2 + a1*a3*b3
    o2 = (-a1*a1 - a3*a3)*b2 + a1*a2*b1 + a2*a3*b3
    o3 = (-a2*a2 - a1*a1)*b3 + a1*a3*b1 + a2*a3*b2
    return o1, o2, o3
end

@inline function evalpoly3(t, C::AbstractVector{Vector3{T}}) where T
    x1 = zero(T)
    x2 = zero(T)
    x3 = zero(T)

    @inbounds for k in length(C):-1:1
        c = C[k]
        x1 = muladd(x1, t, c[1])
        x2 = muladd(x2, t, c[2])
        x3 = muladd(x3, t, c[3])
    end

    return (x1, x2, x3)
end

# VECTOR FUNCTIONS ----------------------------------------------------------------------- #
@inline Base.iszero(a::Vector3) = all(iszero, a)
@inline Base.:+(a::Vector3, b::Vector3) = sum3(a, b)                # Sum
@inline Base.:+(a::Vector3, b::AbstractVector) = sum3(a, b)         # Sum
@inline Base.:+(a::AbstractVector, b::Vector3) = sum3(a, b)         # Sum

@inline Base.:-(v::Vector3) = (-v[1], -v[2], -v[3])                 # Unary minus

@inline Base.:-(a::Vector3, b::Vector3) = subt3(a, b)               # Subtraction
@inline Base.:-(a::Vector3, b::AbstractVector) = subt3(a, b)        # Subtraction
@inline Base.:-(a::AbstractVector, b::Vector3) = subt3(a, b)        # Subtraction

@inline Base.:*(a::Vector3, s::Number) = smult3(s, a)               # Scalar multiplication
@inline Base.:*(s::Number, a::Vector3) = smult3(s, a)               # Scalar multiplication
@inline Base.:/(a::Vector3, s::Number) = smult3(1/s, a)             # Scalar division
@inline Base.zero(::Type{Vector3{T}}) where {T} = (zero(T), zero(T), zero(T))

@inline LinearAlgebra.dot(a::Vector3, b::Vector3) = dot3(a, b)
@inline LinearAlgebra.dot(a::Vector3, b::AbstractVector) = dot3(a, b)
@inline LinearAlgebra.dot(a::AbstractVector, b::Vector3) = dot3(a, b)

@inline LinearAlgebra.cross(a::Vector3, b::Vector3) = cross3(a, b)
@inline LinearAlgebra.cross(a::Vector3, b::AbstractVector) = cross3(a, b)
@inline LinearAlgebra.cross(a::AbstractVector, b::Vector3) = cross3(a, b)
@inline LinearAlgebra.transpose(a::Vector3) = toVector(a)'          # This allocates, but it is sometimes needed
@inline LinearAlgebra.adjoint(a::Vector3) = toVector(a)'            # This allocates, but it is sometimes needed

@inline crossSq(a::AbstractVector, b::AbstractVector) = toVector(crossSq3(a, b))
@inline crossSq(a, b) = crossSq3(a, b)

@inline LinearAlgebra.norm(a::Vector3, p::Real=2) =
    p == 2 ? sqrt(dot(a, a)) :
    p == Inf ? maximum(abs, a) :
    p == 1 ? sum(abs, a) :
    sum(abs(x)^p for x in a)^(1/p)

@inline LinearAlgebra.normalize(a::Vector3) = a / norm(a)

@inline Vector3(v::AbstractVector) = (v[1], v[2], v[3])
@inline toVector(v::Vector3) = [v[1], v[2], v[3]]
@inline toVector!(v::AbstractVector, a::Vector3) = v .= a

# QUATERNION FUNCTIONS ------------------------------------------------------------------- #
@inline Base.iszero(a::Quaternion) = all(iszero, a)
@inline Base.:*(s::Number, a::Quaternion) = (s*a[1], s*a[2], s*a[3], s*a[4])               # Scalar multiplication
@inline Base.:*(a::Quaternion, s::Number) = (s*a[1], s*a[2], s*a[3], s*a[4])               # Scalar multiplication
@inline Base.:/(a::Quaternion, s::Number) = (a[1]/s, a[2]/s, a[3]/s, a[4]/s)               # Scalar division
@inline Base.:*(p::Quaternion, q::Quaternion) = q_multiplyCore(p, q)
@inline Base.:*(v::Vector3{T}, q::Quaternion) where {T} = q_multiplyCore((zero(T), v[1], v[2], v[3]), q)
@inline Base.:*(q::Quaternion, v::Vector3{T}) where {T} = q_multiplyCore(q, (zero(T), v[1], v[2], v[3]))
@inline LinearAlgebra.norm(q::Quaternion) = sqrt(q[1]*q[1] + q[2]*q[2] + q[3]*q[3] + q[4]*q[4])
@inline LinearAlgebra.normalize(q::Quaternion) = q / norm(q)
@inline LinearAlgebra.conj(q::Quaternion) = (q[1], -q[2], -q[3], -q[4])
@inline LinearAlgebra.adjoint(q::Quaternion) = (q[1], -q[2], -q[3], -q[4])

@inline Quaternion(q::AbstractVector) = (q[1], q[2], q[3], q[4])
@inline toVector(q::Quaternion) = [q[1], q[2], q[3], q[4]]
@inline toVector!(v::AbstractVector, q::Quaternion) = v .= q

# MATRIX FUNCTIONS ----------------------------------------------------------------------- #
@inline mrow3(R::AbstractMatrix, i) = (R[i, 1], R[i, 2], R[i, 3])
@inline mcol3(R::AbstractMatrix, i) = (R[1, i], R[2, i], R[3, i])

@inline function mrow3(R::Matrix3x3, i)
    a11, a21, a31, a12, a22, a32, a13, a23, a33 = R
    i == 1 && return (a11, a12, a13)
    i == 2 && return (a21, a22, a23)
    return (a31, a32, a33)
end

@inline function mcol3(R::Matrix3x3, i)
    a11, a21, a31, a12, a22, a32, a13, a23, a33 = R
    i == 1 && return (a11, a21, a31)
    i == 2 && return (a12, a22, a32)
    return (a13, a23, a33)
end

function madd3x3(A, B)
    a11, a21, a31, a12, a22, a32, a13, a23, a33 = A
    b11, b21, b31, b12, b22, b32, b13, b23, b33 = B

    c11 = a11 + b11; c12 = a12 + b12; c13 = a13 + b13
    c21 = a21 + b21; c22 = a22 + b22; c23 = a23 + b23
    c31 = a31 + b31; c32 = a32 + b32; c33 = a33 + b33
    return (c11, c21, c31, c12, c22, c32, c13, c23, c33)
end

function msubt3x3(A, B)
    a11, a21, a31, a12, a22, a32, a13, a23, a33 = A
    b11, b21, b31, b12, b22, b32, b13, b23, b33 = B

    c11 = a11 - b11; c12 = a12 - b12; c13 = a13 - b13
    c21 = a21 - b21; c22 = a22 - b22; c23 = a23 - b23
    c31 = a31 - b31; c32 = a32 - b32; c33 = a33 - b33
    return (c11, c21, c31, c12, c22, c32, c13, c23, c33)
end

function mmult3x3(A, B)
    a11, a21, a31, a12, a22, a32, a13, a23, a33 = A
    b11, b21, b31, b12, b22, b32, b13, b23, b33 = B

    c11 = a11 * b11 + a12 * b21 + a13 * b31
    c12 = a11 * b12 + a12 * b22 + a13 * b32
    c13 = a11 * b13 + a12 * b23 + a13 * b33

    c21 = a21 * b11 + a22 * b21 + a23 * b31
    c22 = a21 * b12 + a22 * b22 + a23 * b32
    c23 = a21 * b13 + a22 * b23 + a23 * b33

    c31 = a31 * b11 + a32 * b21 + a33 * b31
    c32 = a31 * b12 + a32 * b22 + a33 * b32
    c33 = a31 * b13 + a32 * b23 + a33 * b33
    return (c11, c21, c31, c12, c22, c32, c13, c23, c33)
end


@inline Matrix3x3(M::AbstractMatrix) = (M[1], M[2], M[3], M[4], M[5], M[6], M[7], M[8], M[9])
@inline Base.:*(s::Number, M::Matrix3x3) = (s*M[1], s*M[2], s*M[3], s*M[4], s*M[5], s*M[6], s*M[7], s*M[8], s*M[9])               # Scalar multiplication
@inline Base.:*(M::Matrix3x3, s::Number) = s*M                          # Scalar multiplication
@inline Base.:*(M::Matrix3x3, v::Vector3) = mmul3(M, v)                 # Matrix-vector multiplication
@inline Base.:*(M::Matrix3x3, v::AbstractVector) = mmul3(M, v)          # Matrix-vector multiplication
@inline Base.:*(M::AbstractMatrix, v::Vector3) = mmul3(M, v)            # Matrix-vector multiplication
@inline Base.:*(A::Matrix3x3, B::Matrix3x3) = mmult3x3(A, B)            # Matrix-matrix multiplication
@inline Base.:*(A::AbstractMatrix, B::Matrix3x3) = mmult3x3(A, B)       # Matrix-matrix multiplication
@inline Base.:*(A::Matrix3x3, B::AbstractMatrix) = mmult3x3(A, B)       # Matrix-matrix multiplication

@inline Base.:+(A::Matrix3x3, B::Matrix3x3) = madd3x3(A, B)             # Sum
@inline Base.:+(A::Matrix3x3, B::AbstractMatrix) = madd3x3(A, B)        # Sum
@inline Base.:+(A::AbstractMatrix, B::Matrix3x3) = madd3x3(A, B)        # Sum
@inline Base.:+(A::Matrix3x3, B::UniformScaling) = madd3x3(A, B)        # Sum

@inline Base.:-(M::Matrix3x3) = (-M[1], -M[2], -M[3], -M[4], -M[5], -M[6], -M[7], -M[8], -M[9])     # Unary minus
@inline Base.:-(A::Matrix3x3, B::Matrix3x3) = msubt3x3(A, B)            # Subtraction
@inline Base.:-(A::Matrix3x3, B::AbstractMatrix) = msubt3x3(A, B)       # Subtraction
@inline Base.:-(A::AbstractMatrix, B::Matrix3x3) = msubt3x3(A, B)       # Subtraction

@inline function LinearAlgebra.diagm(v::Vector3{T}) where{T}
    x, y, z = v
    o = zero(T)
    return (x, o, o, o, y, o, o, o, z)
end

@inline LinearAlgebra.adjoint(M::Matrix3x3) = (M[1], M[4], M[7], M[2], M[5], M[8], M[3], M[6], M[9])
@inline LinearAlgebra.transpose(M::Matrix3x3) = (M[1], M[4], M[7], M[2], M[5], M[8], M[3], M[6], M[9])
@inline function LinearAlgebra.det(M::Matrix3x3)
    a, d, g, b, e, h, c, f, i = M
    return a * (e * i - f * h) - b * (d * i - g * f) + c * (d * h - e * g)
end

@inline toMatrix(M::Matrix3x3) = [M[1] M[4] M[7]; M[2] M[5] M[8]; M[3] M[6] M[9]]
@inline function toMatrix!(Mout::AbstractMatrix, M::Matrix3x3)
    @inbounds for i in 1:9
        Mout[i] = M[i]
    end
    return Mout
end

getij(M::Matrix3x3, i, j) = M[i + 3(j-1)]


# DISPLAY FUNCTIONS ---------------------------------------------------------------------- #
function Base.show(io::IO, ::MIME"text/plain", v::Vector3{T}) where T
    println(io, "Quats.Vector3{$T}:")
    for i in 1:3
        print(io, " ")
        x = v[i]
        !signbit(x) && print(io, " ")
        show(io, round(x, digits=16))
        i < 3 && println(io)
    end
end

function Base.show(io::IO, ::MIME"text/plain", q::Quaternion{T}) where T
    println(io, "Quats.Quaternion{$T}:")
    for i in 1:4
        print(io, " ")
        x = q[i]
        !signbit(x) && print(io, " ")
        show(io, round(x, digits=16))
        i < 4 && println(io)
    end
end

function Base.show(io::IO, ::MIME"text/plain", M::Matrix3x3{T}) where T
    println(io, "Quats.Matrix3x3{$T}:")
    for i in 1:3
        print(io, " ")
        for j in 1:3
            j > 1 && print(io, "  ")
            x = M[i + 3(j-1)]
            !signbit(x) && print(io, " ")
            show(io, round(x, digits=6))
        end
        i < 3 && println(io)
    end
end
