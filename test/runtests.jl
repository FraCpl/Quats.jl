using Quats
using LinearAlgebra
using Test

function TEST_rots()
    N = 100
    v_A = [randn(3) for _ in 1:N]
    vq_A = [[0.0; v_A[i]] for i in 1:N]
    R_BA = [dcm_random() for _ in 1:N]
    q_BA = q_fromDcm.(R_BA)
    q_AB = q_transpose.(q_BA)
    v1_B = [R_BA[i]*v_A[i] for i in 1:N]
    v2_B = q_transformVector.(q_BA, v_A)
    vq3_B = q_multiply.(q_BA, q_multiply.(vq_A, q_AB))
    v3_B = [vq3_B[i][2:4] for i in 1:N]

    axis = normalize(randn(3))
    angle = 2π*rand()

    R1 = q_toDcm(q_fromAxisAngle(axis, angle))
    R2 = dcm_fromAxisAngle(axis, angle)

    ε = maximum([
        maximum(norm.(v1_B - v2_B));
        maximum(norm.(v1_B - v3_B));
        maximum(norm.(q_toDcm.(q_BA) - R_BA));
        maximum(norm.(R1 - R2))
    ])
    return ε
end

function TEST_qMultiplyn()
    q = [q_random() for _ in 1:4]
    qTrue = q_multiply(q[1], q_multiply(q[2], q_multiply(q[3], q[4])))
    qOut = q_random()

    err1 = q_multiplyn(q...) - qTrue
    q_multiplyn!(qOut, q...)
    err2 = qOut - qTrue
    return max(maximum(abs, err1), maximum(abs, err2))
end

function TEST_qMultiply()
    p = q_random()
    q = q_random()

    ps = p[1]; pv = p[2:4]
    qs = q[1]; qv = q[2:4]

    pq = q_multiply(p, q)
    pqTrue = [ps*qs - dot(pv, qv); ps*qv + qs*pv + cross(pv, qv)]
    return norm(pq - pqTrue)
end

function TEST_qFromDcm()
    qTrue = q_random()
    Rtrue = q_toDcm(qTrue)
    q = q_fromDcm(Rtrue)
    return norm(q - sign(q[1]*qTrue[1])*qTrue)
end

function TEST_qTransformVector()
    q = q_random()
    v = randn(3)
    u = q_transformVector(q, v)
    uTrue = q_multiplyn(q, [0.0; v], q_transpose(q))[2:4]
    return norm(u - uTrue)
end

function TEST_qDerivative()
    q = q_random()
    w = randn(3)
    qDot = q_derivative(q, w)
    qDotTrue = q_multiply(q, [0; w]) ./ 2
    return norm(qDot - qDotTrue)
end

#=
function TEST_kin()
    ω0 = 10.0*π/180
    Δt = 23.6
    ε = zeros(6)

    for i in 1:6
        ω = zeros(3)
        if i < 4
            ω[i] = ω0
        else
            ω[i-3] = -ω0
        end
        sol = solve(ODEProblem((q, p, t) -> q_derivative(q, ω),[1.0; 0.0; 0.0; 0.0], (0.0, Δt)); reltol=1e-12, abstol=1e-12)
        qf = q_fromAxisAngle(normalize(ω), ω0*Δt)
        ε[i] = maximum(abs.(q_multiply(sol.u[end],q_transpose(qf))[2:4]))
    end

    return maximum(ε)
end
=#

function TEST_euler123()
    seq = :xyz
    eul = [[-π + 2π*rand(), -π/2 + π*rand(), -π + 2π*rand()] for _ in 1:10_000]
    eult = [collect(q_toEuler(q_fromEuler(ek..., seq), seq)) for ek in eul]
    return maximum(norm.(eul - eult))
end

function TEST_euler321()
    seq = :zyx
    eul = [[-π + 2π*rand(), -π/2 + π*rand(), -π + 2π*rand()] for _ in 1:10_000]
    eult = [collect(q_toEuler(q_fromEuler(ek..., seq), seq)) for ek in eul]
    return maximum([norm(eul[k] - eult[k]) for k in eachindex(eul)])
end

function TEST_euler132()
    seq = :xzy
    eul = [[-π + 2π*rand(), -π/2 + π*rand(), -π + 2π*rand()] for _ in 1:10_000]
    eult = [collect(q_toEuler(q_fromEuler(ek..., seq), seq)) for ek in eul]
    return maximum([norm(eul[k] - eult[k]) for k in eachindex(eul)])
end

function TEST_euler213()
    seq = :yxz
    eul = [[-π + 2π*rand(), -π/2 + π*rand(), -π + 2π*rand()] for _ in 1:10_000]
    eult = [collect(q_toEuler(q_fromEuler(ek..., seq), seq)) for ek in eul]
    return maximum([norm(eul[k] - eult[k]) for k in eachindex(eul)])
end

# function TEST_euler231()
#     seq = :yzx
#     eul = [[-π + 2π*rand(), -π + 2π*rand(), -π/2 + π*rand()] for _ in 1:10_000]
#     eult = [collect(q_toEuler(q_fromEuler(ek..., seq), seq)) for ek in eul]
#     return maximum([norm(eul[k] - eult[k]) for k in eachindex(eul)])
# end

function TEST_cross()
    a = randn(3)
    b = randn(3)
    x = zeros(3)
    cross!(x, a, b)
    err1 = norm(x - cross(a, b))
    crossSq!(x, a, b)
    err2 = norm(x - cross(a, cross(a, b)))
    return maximum([err1; err2])
end

function TEST_qtoxyz()
    q = q_random()
    ex = norm(q_tox(q) - q_transformVector(q, [1.0; 0.0; 0.0]))
    ey = norm(q_toy(q) - q_transformVector(q, [0.0; 1.0; 0.0]))
    ez = norm(q_toz(q) - q_transformVector(q, [0.0; 0.0; 1.0]))
    return maximum([ex; ey; ez])
end

function TEST_qToDcm()
    q = q_random()
    R = q_toDcm(q)
    qvx = crossMat(q[2:4])
    Rtrue = I + 2*q[1].*qvx + 2*qvx*qvx
    return norm(R*Rtrue' - I)
end

@testset "Quats.jl" begin
    ERR_TOL = 1e-10
    @test TEST_rots() < ERR_TOL
    @test TEST_qMultiply() < ERR_TOL
    @test TEST_qMultiplyn() < ERR_TOL
    @test TEST_euler123() < ERR_TOL
    @test TEST_euler321() < ERR_TOL
    @test TEST_euler132() < ERR_TOL
    @test TEST_euler213() < ERR_TOL
    # @test TEST_euler231() < ERR_TOL
    @test TEST_cross() < ERR_TOL
    #@test TEST_kin() < ERR_TOL
    @test TEST_qtoxyz() < ERR_TOL
    @test TEST_qToDcm() < ERR_TOL
    @test TEST_qFromDcm() < ERR_TOL
    @test TEST_qTransformVector() < ERR_TOL
    @test TEST_qDerivative() < ERR_TOL
end
