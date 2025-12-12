```@meta
CurrentModule = Quats
```

# Quats.jl

This package provides a bare-bones representation of attitude transformations through quaternions, rotation
vectors, and transformation matrices (a.k.a. direction cosine matrices -- DCM). 

# Installation
Just type:
```julia
using Pkg
Pkg.add(url="https://github.com/FraCpl/Quats.jl")
```

# Notation and conventions
Given a 3-dimensional vector ``v``, a reference frame ``A`` and a reference frame ``B``, we indicate with 
``v^A`` the projection of ``v`` into ``A``, and with ``v^B`` its projection in ``B``, so that the two 
projections can be related through the transformation matrix ``R_{AB}`` as follows
```math
v^A = R_{AB} v^B
```
The same transformation can be represented using the quaternion ``q_{AB}``
```math
\begin{bmatrix} 0\\ v^A\end{bmatrix} = q_{AB} ⊗ \begin{bmatrix} 0\\ v^B\end{bmatrix} ⊗ q_{AB}^*
```

In this package, quaternions are represented as 4-dimensional vectors, using the basic vector format
natively provided by Julia. The scalar component of the quaternion corresponds to its first element, 
i.e., ```q[1]```. We use a right-handed passive convention for representing attitude 
transformations through quaternions, and use ``\theta_{AB}`` to indicate the angle by which the frame 
``A`` needs to be actively rotated to make it coincide with frame ``B`` (i.e., the rotation from ``A``
to ``B``).

!!! note
    Other (and weirder) quaternions notations exist in the litterature, which are obviously **not**
    covered by this package.

# API

### Quats
```@autodocs
Modules = [Quats]
Order   = [:function, :type]
Pages   = ["q.jl"]
```
### Transformation matrices
```@autodocs
Modules = [Quats]
Order   = [:function, :type]
Pages   = ["dcm.jl"]
```
### Rotation vectors
```@autodocs
Modules = [Quats]
Order   = [:function, :type]
Pages   = ["rv.jl"]
```
### Cross product functions
```@autodocs
Modules = [Quats]
Order   = [:function, :type]
Pages   = ["cross.jl"]
```
## Index
```@index
```