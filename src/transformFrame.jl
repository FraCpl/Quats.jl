# Provides pos_A, vel_A. Outputs are provided as Vector3
function transformFrame(pos_B, vel_B, q_AB, ωAB_B)
    pos_A = Quats.q_transformVectorCore(q_AB, pos_B)
    vel_A = Quats.q_transformVectorCore(q_AB, vel_B + Quats.cross3(ωAB_B, pos_B))
    return pos_A, vel_A
end

# Proivdes posAX_A, velAX_A. Outputs are provided as Vector3
function transformFrame(posBX_B, velBX_B, q_AB, ωAB_B, posAB_A, velAB_A)
    posBX_A, velBX_A = transformFrame(posBX_B, velBX_B, q_AB, ωAB_B)
    return posAB_A + posBX_A, velAB_A + velBX_A
end
