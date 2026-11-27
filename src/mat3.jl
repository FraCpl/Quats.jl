function crossMat3core(v)
    x, y, z = v
    o = zero(x)
    return (o, z, -y, -z, o, x, y, -x, o)
end

function crossMatSq3core(v)
    x, y, z = v
    x2 = x*x; y2 = y*y; z2 = z*z
    return (-y2 - z2, x*y, x*z, x*y, -x2 - z2, y*z, x*z, y*z, -x2 - y2)
end

"""
    [v×] = crossMat(v)

Compute the cross product matrix of a vector ```v```.
"""
@inline crossMat(v::AbstractArray) = toMatrix(crossMat3core(v))
@inline crossMat(v::Vec3) = crossMat3core(v)
@inline crossMat!(R::AbstractMatrix, v) = toMatrix!(R, crossMat3core(v))

"""
    [v×]² = crossMatSq(x, y, z)

Compute the squared cross product matrix of a vector ```v``` given its scalar components.
"""
@inline crossMatSq(v::AbstractArray) = toMatrix(crossMatSq3core(v))
@inline crossMatSq(v::Vec3) = crossMatSq3core(v)
@inline crossMatSq!(R::AbstractMatrix, v) = toMatrix!(R, crossMatSq3core(v))


# MATRIX FUNCTIONS ----------------------------------------------------------------------- #
@inline mrow3(R::AbstractMatrix, i) = (R[i, 1], R[i, 2], R[i, 3])
@inline mcol3(R::AbstractMatrix, i) = (R[1, i], R[2, i], R[3, i])

@inline function mrow3(R::Quats.Mat3, i)
    a11, a21, a31, a12, a22, a32, a13, a23, a33 = R
    i == 1 && return a11, a12, a13
    i == 2 && return a21, a22, a23
    return a31, a32, a33
end

@inline function mcol3(R::Quats.Mat3, i)
    a11, a21, a31, a12, a22, a32, a13, a23, a33 = R
    i == 1 && return a11, a21, a31
    i == 2 && return a12, a22, a32
    return a13, a23, a33
end

# A + B
function madd3x3(A, B)
    a11, a21, a31, a12, a22, a32, a13, a23, a33 = A
    b11, b21, b31, b12, b22, b32, b13, b23, b33 = B

    c11 = a11 + b11; c12 = a12 + b12; c13 = a13 + b13
    c21 = a21 + b21; c22 = a22 + b22; c23 = a23 + b23
    c31 = a31 + b31; c32 = a32 + b32; c33 = a33 + b33
    return c11, c21, c31, c12, c22, c32, c13, c23, c33
end

# A - B
function msubt3x3(A, B)
    a11, a21, a31, a12, a22, a32, a13, a23, a33 = A
    b11, b21, b31, b12, b22, b32, b13, b23, b33 = B

    c11 = a11 - b11; c12 = a12 - b12; c13 = a13 - b13
    c21 = a21 - b21; c22 = a22 - b22; c23 = a23 - b23
    c31 = a31 - b31; c32 = a32 - b32; c33 = a33 - b33
    return c11, c21, c31, c12, c22, c32, c13, c23, c33
end

# A * B
function mmul3x3(A, B)
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
    return c11, c21, c31, c12, c22, c32, c13, c23, c33
end

# Constructor from AbstractMatrix
@inline Quats.Mat3(M::AbstractMatrix) = M[1], M[2], M[3], M[4], M[5], M[6], M[7], M[8], M[9]

# Scalar multiplication
@inline Base.:*(s::Number, M::Quats.Mat3) = s*M[1], s*M[2], s*M[3], s*M[4], s*M[5], s*M[6], s*M[7], s*M[8], s*M[9]
@inline Base.:*(M::Quats.Mat3, s::Number) = s*M

# Matrix-vector multiplication
@inline Base.:*(M::Quats.Mat3, v::Vec3) = mmul3(M, v)
@inline Base.:*(M::Quats.Mat3, v::AbstractVector) = mmul3(M, v)
@inline Base.:*(M::AbstractMatrix, v::Vec3) = mmul3(M, v)

# Transpose vector-matrix multiplication
@inline Base.:*(v::Vec3, M::Quats.Mat3) = mmul3t(v, M)
@inline Base.:*(v::AbstractVector, M::Quats.Mat3) = mmul3t(v, M)
@inline Base.:*(v::Vec3, M::AbstractMatrix) = mmul3t(v, M)

# Matrix-matrix multiplication
@inline Base.:*(A::Quats.Mat3, B::Quats.Mat3) = mmul3x3(A, B)
@inline Base.:*(A::AbstractMatrix, B::Quats.Mat3) = mmul3x3(A, B)
@inline Base.:*(A::Quats.Mat3, B::AbstractMatrix) = mmul3x3(A, B)

# Sum
@inline Base.:+(A::Quats.Mat3, B::Quats.Mat3) = madd3x3(A, B)
@inline Base.:+(A::Quats.Mat3, B::AbstractMatrix) = madd3x3(A, B)
@inline Base.:+(A::AbstractMatrix, B::Quats.Mat3) = madd3x3(A, B)

# Unary minus
@inline Base.:-(M::Quats.Mat3) = -M[1], -M[2], -M[3], -M[4], -M[5], -M[6], -M[7], -M[8], -M[9]

# Subtraction
@inline Base.:-(A::Quats.Mat3, B::Quats.Mat3) = msubt3x3(A, B)
@inline Base.:-(A::Quats.Mat3, B::AbstractMatrix) = msubt3x3(A, B)
@inline Base.:-(A::AbstractMatrix, B::Quats.Mat3) = msubt3x3(A, B)

@inline function LinearAlgebra.diagm(v::Vec3{T}) where{T}
    x, y, z = v
    o = zero(T)
    return x, o, o, o, y, o, o, o, z
end

@inline LinearAlgebra.adjoint(M::Quats.Mat3) = M[1], M[4], M[7], M[2], M[5], M[8], M[3], M[6], M[9]
@inline LinearAlgebra.transpose(M::Quats.Mat3) = M[1], M[4], M[7], M[2], M[5], M[8], M[3], M[6], M[9]
@inline function LinearAlgebra.det(M::Quats.Mat3)
    a, d, g, b, e, h, c, f, i = M
    return a * (e * i - f * h) - b * (d * i - g * f) + c * (d * h - e * g)
end

@inline toMatrix(M::Quats.Mat3) = [M[1] M[4] M[7]; M[2] M[5] M[8]; M[3] M[6] M[9]]

@inline function Base.copyto!(Mout::AbstractMatrix, M::Quats.Mat3)
    @inbounds for i in 1:9
        Mout[i] = M[i]
    end
    return Mout
end

getij(M::Quats.Mat3, i, j) = M[i + 3(j-1)]

# DISPLAY FUNCTIONS ---------------------------------------------------------------------- #
function Base.show(io::IO, ::MIME"text/plain", M::Quats.Mat3{T}) where T
    println(io, "Quats.Quats.Mat3{$T}:")
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
