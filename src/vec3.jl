# CORE COMPONENTS ======================================================================== #
@inline zeros3() = (0.0, 0.0, 0.0)
@inline randn3() = (randn(), randn(), randn())
@inline randn3(rng) = (randn(rng), randn(rng), randn(rng))
@inline rand3() = (rand(), rand(), rand())
@inline rand3(rng) = (rand(rng), rand(rng), rand(rng))

# Basic math on vectors
@inline sum3(a, b) = a[1]+b[1], a[2]+b[2], a[3]+b[3]
@inline subt3(a, b) = a[1]-b[1], a[2]-b[2], a[3]-b[3]
@inline smult3(s, a) = s*a[1], s*a[2], s*a[3]
@inline dot3(a, b) = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
@inline cross3(a, b) = a[2]*b[3] - a[3]*b[2], a[3]*b[1] - a[1]*b[3], a[1]*b[2] - a[2]*b[1]
@inline vmul3core(a, b) = a[1]*b[1], a[2]*b[1], a[3]*b[1], a[1]*b[2], a[2]*b[2], a[3]*b[2], a[1]*b[3], a[2]*b[3], a[3]*b[3]   # Dyadic product of two vectors, i.e., a * b'

# R * v where R is a 3x3 matrix, and v is a 3x1 vector
@inline function mmul3(R, v)
    r11, r21, r31, r12, r22, r32, r13, r23, r33 = R
    x, y, z = v
    ox = r11 * x + r12 * y + r13 * z
    oy = r21 * x + r22 * y + r23 * z
    oz = r31 * x + r32 * y + r33 * z
    return ox, oy, oz
end
@inline mmul3!(out, R, v) = out .= mmul3(R, v)

# v' * R, where R is a 3x3 matrix, and v is a 3x1 vector
@inline function mmul3t(v, R)
    r11, r21, r31, r12, r22, r32, r13, r23, r33 = R
    x, y, z = v
    ox = x * r11 + y * r21 + z * r31
    oy = x * r12 + y * r22 + z * r32
    oz = x * r13 + y * r23 + z * r33
    return ox, oy, oz
end

# Double cross product: a × (a × b)
function crossSq3(a, b)
    a1, a2, a3 = a
    b1, b2, b3 = b
    o1 = (-a2*a2 - a3*a3)*b1 + a1*a2*b2 + a1*a3*b3
    o2 = (-a1*a1 - a3*a3)*b2 + a1*a2*b1 + a2*a3*b3
    o3 = (-a2*a2 - a1*a1)*b3 + a1*a3*b1 + a2*a3*b2
    return o1, o2, o3
end

# Polynomial interpolation (slightly faster than 'evalpoly')
@inline function evalpoly3(t, C::AbstractVector{Quats.Vec3{T}}) where T
    x = y = z = zero(T)
    @inbounds for k in length(C):-1:1
        cx, cy, cz = C[k]
        x = muladd(x, t, cx)
        y = muladd(y, t, cy)
        z = muladd(z, t, cz)
    end
    return x, y, z
end

# API ==================================================================================== #
# Constructor
@inline Quats.Vec3(v::AbstractVector) = v[1], v[2], v[3]

# Zeros
@inline Base.zero(::Type{Quats.Vec3{T}}) where {T} = zero(T), zero(T), zero(T)
@inline Base.iszero(a::Quats.Vec3) = all(iszero, a)

# From/to Vector
@inline toVector(v::Quats.Vec3) = [v[1], v[2], v[3]]
@inline Base.copyto!(v::AbstractVector, a::Quats.Vec3) = v .= a

# Sum
@inline Base.:+(a::Quats.Vec3, b::Quats.Vec3) = sum3(a, b)
@inline Base.:+(a::Quats.Vec3, b::AbstractVector) = sum3(a, b)
@inline Base.:+(a::AbstractVector, b::Quats.Vec3) = sum3(a, b)

# Unary minus
@inline Base.:-(v::Quats.Vec3) = -v[1], -v[2], -v[3]

# Subtraction
@inline Base.:-(a::Quats.Vec3, b::Quats.Vec3) = subt3(a, b)
@inline Base.:-(a::Quats.Vec3, b::AbstractVector) = subt3(a, b)
@inline Base.:-(a::AbstractVector, b::Quats.Vec3) = subt3(a, b)

# Scalar operations
@inline Base.:*(a::Quats.Vec3, s::Number) = smult3(s, a)
@inline Base.:*(s::Number, a::Quats.Vec3) = smult3(s, a)
@inline Base.:/(a::Quats.Vec3, s::Number) = smult3(1/s, a)

# Dot product
@inline LinearAlgebra.dot(a::Quats.Vec3, b::Quats.Vec3) = dot3(a, b)
@inline LinearAlgebra.dot(a::Quats.Vec3, b::AbstractVector) = dot3(a, b)
@inline LinearAlgebra.dot(a::AbstractVector, b::Quats.Vec3) = dot3(a, b)

# Cross product
@inline LinearAlgebra.cross(a::Quats.Vec3, b::Quats.Vec3) = cross3(a, b)
@inline LinearAlgebra.cross(a::Quats.Vec3, b::AbstractVector) = cross3(a, b)
@inline LinearAlgebra.cross(a::AbstractVector, b::Quats.Vec3) = cross3(a, b)
@inline cross!(a, b, c) = a .= cross3(b, c)
@inline addCross!(a, b, c) = a .+= cross3(b, c)

@inline crossSq(a::AbstractVector, b::AbstractVector) = toVector(crossSq3(a, b))
@inline crossSq(a, b) = crossSq3(a, b)
@inline crossSq!(a, b, c) = a .= crossSq3(b, c)
@inline addCrossSq!(a, b, c) = a .+= crossSq3(b, c)

# Transpose
# This allocates, but it is sometimes needed
@inline LinearAlgebra.transpose(a::Quats.Vec3) = toVector(a)'
@inline LinearAlgebra.adjoint(a::Quats.Vec3) = toVector(a)'

# Norm & Normalization
@inline LinearAlgebra.norm(a::Quats.Vec3) = hypot(a[1], a[2], a[3])
@inline LinearAlgebra.normalize(a::Quats.Vec3) = a / norm(a)

# Dyadic product of two vectors, i.e., returning a 3x3 matrix given by a * b'
@inline vmul3(a::AbstractVector, b::AbstractVector=a) = toMatrix(vmul3core(a, b))
@inline vmul3(a, b=a) = vmul3core(a, b)
@inline vmul3!(out, a, b=a) = out .= vmul3core(a, b)

# DISPLAY FUNCTIONS ---------------------------------------------------------------------- #
function Base.show(io::IO, ::MIME"text/plain", v::Vec3{T}) where T
    println(io, "Quats.Vec3{$T}:")
    for i in 1:3
        print(io, " ")
        x = v[i]
        !signbit(x) && print(io, " ")
        show(io, round(x, digits=16))
        i < 3 && println(io)
    end
end
