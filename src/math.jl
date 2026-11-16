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
@inline mrow3(R, i) = (R[i, 1], R[i, 2], R[i, 3])
@inline mcol3(R, i) = (R[1, i], R[2, i], R[3, i])

function crossSq3(a, b)
    a1, a2, a3 = a
    b1, b2, b3 = b
    o1 = (-a2*a2 - a3*a3)*b1 + a1*a2*b2 + a1*a3*b3
    o2 = (-a1*a1 - a3*a3)*b2 + a1*a2*b1 + a2*a3*b3
    o3 = (-a2*a2 - a1*a1)*b3 + a1*a3*b1 + a2*a3*b2
    return o1, o2, o3
end

const Vector3{T} = NTuple{3, T}
@inline Base.:+(a::Vector3, b::Vector3) = sum3(a, b)         # Sum
@inline Base.:+(a::Vector3, b::AbstractVector) = sum3(a, b)  # Sum
@inline Base.:+(a::AbstractVector, b::Vector3) = sum3(a, b)  # Sum

@inline Base.:-(v::Vector3) = (-v[1], -v[2], -v[3])                   # Unary minus

@inline Base.:-(a::Vector3, b::Vector3) = subt3(a, b)         # Subtraction
@inline Base.:-(a::Vector3, b::AbstractVector) = subt3(a, b)  # Subtraction
@inline Base.:-(a::AbstractVector, b::Vector3) = subt3(a, b)  # Subtraction

@inline Base.:*(a::Vector3, s::Number) = smult3(s, a)              # Scalar multiplication
@inline Base.:*(s::Number, a::Vector3) = smult3(s, a)              # Scalar multiplication
@inline Base.:/(a::Vector3, s::Number) = smult3(1/s, a)            # Scalar division
@inline Base.zero(::Type{Vector3{T}}) where {T} = (zero(T), zero(T), zero(T))

@inline LinearAlgebra.dot(a::Vector3, b::Vector3) = dot3(a, b)
@inline LinearAlgebra.dot(a::Vector3, b::AbstractVector) = dot3(a, b)
@inline LinearAlgebra.dot(a::AbstractVector, b::Vector3) = dot3(a, b)

@inline LinearAlgebra.cross(a::Vector3, b::Vector3) = cross3(a, b)
@inline LinearAlgebra.cross(a::Vector3, b::AbstractVector) = cross3(a, b)
@inline LinearAlgebra.cross(a::AbstractVector, b::Vector3) = cross3(a, b)

@inline crossSq(a::Vector3, b::Vector3) = crossSq3(a, b)
@inline crossSq(a::Vector3, b::AbstractVector) = crossSq3(a, b)
@inline crossSq(a::AbstractVector, b::Vector3) = crossSq3(a, b)

@inline LinearAlgebra.norm(a::Vector3, p::Real=2) =
    p == 2 ? sqrt(dot(a, a)) :
    p == Inf ? maximum(abs, a) :
    p == 1 ? sum(abs, a) :
    sum(abs(x)^p for x in a)^(1/p)
@inline LinearAlgebra.normalize(a::Vector3) = a / norm(a)

@inline Vector3(v::AbstractVector) = (v[1], v[2], v[3])
@inline toVector(v::Vector3) = collect(v)
@inline toVector!(v::AbstractVector, a::Vector3) = v .= a

const Quaternion{T} = NTuple{4, T}
@inline Base.:*(a::Quaternion, s::Number) = (s*a[1], s*a[2], s*a[3], s*a[4])               # Scalar multiplication
@inline Base.:*(s::Number, a::Quaternion) = (s*a[1], s*a[2], s*a[3], s*a[4])               # Scalar multiplication
@inline Base.:/(a::Quaternion, s::Number) = (a[1]/s, a[2]/s, a[3]/s, a[4]/s)               # Scalar division
@inline Base.:*(p::Quaternion, q::Quaternion) = Quats.q_multiplyCore(p, q)
@inline Base.:*(v::Vector3{T}, q::Quaternion) where {T} = Quats.q_multiplyCore((zero(T), v[1], v[2], v[3]), q)
@inline Base.:*(q::Quaternion, v::Vector3{T}) where {T} = Quats.q_multiplyCore(q, (zero(T), v[1], v[2], v[3]))
@inline LinearAlgebra.norm(q::Quaternion) = sqrt(q[1]*q[1] + q[2]*q[2] + q[3]*q[3] + q[4]*q[4])
@inline LinearAlgebra.normalize(q::Quaternion) = q / norm(q)
@inline LinearAlgebra.conj(q::Quaternion) = (q[1], -q[2], -q[3], -q[4])
@inline LinearAlgebra.adjoint(q::Quaternion) = (q[1], -q[2], -q[3], -q[4])

@inline Quaternion(q::AbstractVector) = (q[1], q[2], q[3], q[4])
@inline toVector(q::Quaternion) = collect(q)
@inline toVector!(v::AbstractVector, q::Quaternion) = v .= q
