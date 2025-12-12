"""
    [v×] = crossMat(v)

Compute the cross product matrix of a vector ```v```.
"""
@inline crossMat(v) = crossMat(v[1], v[2], v[3])

"""
    [v×] = crossMat(x, y, z)

Compute the cross product matrix of a vector ```v``` given its scalar components.
"""
@inline function crossMat(x, y, z)
    R = Matrix{typeof(x)}(undef, 3, 3)
    crossMat!(R, x, y, z)
    return R
end

@inline function crossMat!(R, v)
    crossMat!(R, v[1], v[2], v[3])
    return nothing
end

@inline function crossMat!(R, x, y, z)
    @inbounds for i in 1:3
        ;
        R[i, i] = 0.0;
    end
    R[1, 2] = -z
    R[1, 3] = +y
    R[2, 1] = +z
    R[2, 3] = -x
    R[3, 1] = -y
    R[3, 2] = +x
    return nothing
end

# @inline function crossMat(v::SVector{3, T}) where T
#     return crossMatStatic(v[1], v[2], v[3])
# end

# @inline function crossMatStatic(x, y, z)
#     return @SMatrix [
#         0 -z y; z 0 -x; -y x 0
#     ]
# end

"""
    [v×]² = crossMatSq(x, y, z)

Compute the squared cross product matrix of a vector ```v``` given its scalar components.
"""
@inline function crossMatSq(x, y, z)
    R = Matrix{typeof(x)}(undef, 3, 3)
    crossMatSq!(R, x, y, z)
    return R
end

"""
    [v×]² = crossMatSq(v)

Compute the squared cross product matrix of a vector ```v```.
"""
@inline function crossMatSq(v)
    R = Matrix{eltype(v)}(undef, 3, 3)
    crossMatSq!(R, v[1], v[2], v[3])
    return R
end

@inline function crossMatSq!(R, v)
    crossMatSq!(R, v[1], v[2], v[3])
    return nothing
end

@inline function crossMatSq!(R, x, y, z)
    R[1, 1] = -y*y - z*z
    R[1, 2] = x*y
    R[1, 3] = x*z
    R[2, 2] = -x*x - z*z
    R[2, 3] = y*z
    R[3, 3] = -x*x - y*y
    R[2, 1] = R[1, 2]
    R[3, 1] = R[1, 3]
    R[3, 2] = R[2, 3]
    return nothing
end

# @inline function crossMatSq(v::SVector{3, T}) where T
#     return crossMatSqStatic(v[1], v[2], v[3])
# end

# @inline function crossMatSqStatic(x, y, z)
#     xy = x*y; xz = x*z; yz = y*z
#     xx = x*x; yy = y*y; zz = z*z
#     return @SMatrix [
#         (-yy - zz) xy xz;
#         xy (-xx - zz) yz;
#         xz yz (-xx - yy);
#     ]
# end

@inline crossMatInv(M::Matrix) = [-M[2, 3]; M[1, 3]; -M[1, 2]]
# @inline function crossMatInv(M::SMatrix{3, 3, T}) where T
#     return @SVector [-M[2, 3]; M[1, 3]; -M[1, 2]]
# end

"""
    cross!(a, b, c)

Compute a = b × c.
"""
@inline function cross!(a, b, c)
    a[1] = -b[3]*c[2] + b[2]*c[3]
    a[2] = +b[3]*c[1] - b[1]*c[3]
    a[3] = -b[2]*c[1] + b[1]*c[2]
    return nothing
end

"""
    addCross!(a, b, c)

Compute a += b × c.
"""
@inline function addCross!(a, b, c)
    a[1] += -b[3]*c[2] + b[2]*c[3]
    a[2] += +b[3]*c[1] - b[1]*c[3]
    a[3] += -b[2]*c[1] + b[1]*c[2]
    return nothing
end

"""
    crossSq!(a, b, c)

Compute a = b × (b × c).
"""
@inline function crossSq!(a, b, c)
    b1, b2, b3 = b
    c1, c2, c3 = c

    a[1] = (-b2*b2 - b3*b3)*c1 + b1*b2*c2 + b1*b3*c3
    a[2] = (-b1*b1 - b3*b3)*c2 + b1*b2*c1 + b2*b3*c3
    a[3] = (-b2*b2 - b1*b1)*c3 + b1*b3*c1 + b2*b3*c2
    return nothing
end

"""
    addCrossSq!(a, b, c)

Compute a += b × (b × c).
"""
@inline function addCrossSq!(a, b, c)
    a[1] += (-b2*b2 - b3*b3)*c1 + b1*b2*c2 + b1*b3*c3
    a[2] += (-b1*b1 - b3*b3)*c2 + b1*b2*c1 + b2*b3*c3
    a[3] += (-b2*b2 - b1*b1)*c3 + b1*b3*c1 + b2*b3*c2
    return nothing
end
