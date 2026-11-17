function crossMatCore(v)
    x, y, z = v
    o = zero(x)
    return (o, z, -y, -z, o, x, y, -x, o)
end

function crossMatSqCore(v)
    x, y, z = v
    x2 = x*x; y2 = y*y; z2 = z*z
    return (-y2 - z2, x*y, x*z, x*y, -x2 - z2, y*z, x*z, y*z, -x2 - y2)
end

"""
    [v×] = crossMat(v)

Compute the cross product matrix of a vector ```v```.
"""
@inline crossMat(v::AbstractArray) = toMatrix(crossMatCore(v))
@inline crossMat(v::Vector3) = crossMatCore(v)
@inline crossMat!(R::AbstractMatrix, v) = toMatrix!(R, crossMatCore(v))

"""
    [v×]² = crossMatSq(x, y, z)

Compute the squared cross product matrix of a vector ```v``` given its scalar components.
"""
@inline crossMatSq(v::AbstractArray) = toMatrix(crossMatSqCore(v))
@inline crossMatSq(v::Vector3) = crossMatSqCore(v)
@inline crossMatSq!(R::AbstractMatrix) = toMatrix!(R, crossMatSqCore(v))

"""
    cross!(a, b, c)

Compute a = b × c.
"""
@inline cross!(a, b, c) = a .= cross3(b, c)

"""
    addCross!(a, b, c)

Compute a = a + b × c.
"""
@inline addCross!(a, b, c) = a .+= cross3(b, c)

"""
    crossSq!(a, b, c)

Compute a = b × (b × c).
"""
@inline crossSq!(a, b, c) = a .= crossSq3(b, c)

"""
    addCrossSq!(a, b, c)

Compute a = a + b × (b × c).
"""
@inline addCrossSq!(a, b, c) = a .+= crossSq3(b, c)
