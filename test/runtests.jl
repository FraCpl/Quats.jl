using Quats, LinearAlgebra
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

function TEST_mult()
    q = [q_random() for _ in 1:4]
    ε = q_multiplyn(q[1], q[2], q[3], q[4]) - q_multiply(q[1], q_multiply(q[2], q_multiply(q[3], q[4])))
    return maximum(abs.(ε))
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
    eul = [[-π + 2π*rand(); -π/2 + π*rand(); -π + 2π*rand()] for _ in 1:10_000]
    seq = [1, 2, 3]
    eult = q_toEuler.(q_fromEuler.(eul, Ref(seq)), Ref(seq))
    return maximum(norm.(eul - eult))
end

function TEST_euler321()
    eul = [[-π + 2π*rand(); -π/2 + π*rand(); -π + 2π*rand()] for _ in 1:10_000]
    seq = [3, 2, 1]
    eult = q_toEuler.(q_fromEuler.(eul, Ref(seq)), Ref(seq))
    return maximum(norm.(eul - eult))
end

function TEST_cross()
    a = randn(3)
    b = randn(3)
    x = zeros(3)
    @time cross!(x, a, b)
    err1 = norm(x - cross(a, b))
    @time crossSq!(x, a, b)
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

@testset "Quaternions.jl" begin
    ERR_TOL = 1e-10
    @test TEST_rots() < ERR_TOL
    @test TEST_mult() < ERR_TOL
    @test TEST_euler123() < ERR_TOL
    @test TEST_euler321() < ERR_TOL
    @test TEST_cross() < ERR_TOL
    #@test TEST_kin() < ERR_TOL
    @test TEST_qtoxyz() < ERR_TOL
end
