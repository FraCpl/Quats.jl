
# This function determines the convention used by a quaternion multiplication function in
# terms of: 1) Scalar-first or scalar-last, 2) Right-handed or left-handed algebra
# The multiplication function shall provide q*p = fmult(q, p)
function q_testConvention(fmult)
    q = [1.0; 2.0; 3.0; 4.0]
    p = [-4.0; 2.0; -3.0; 1.0]
    q_testConvention(fmult(q, p), q, p)
end

function q_testConvention(qp, q, p)

    # out = q*p
    function q_multiplyUniversal(q, p; scalarFirst=true, right=true)
        qq = copy(q)
        pp = copy(p)
        if !scalarFirst
            qq = [q[4]; q[1:3]]
            pp = [p[4]; p[1:3]]
        end
        qqpp = right ? qp_right(qq, pp) : qp_left(qq, pp)
        scalarFirst && qqpp;
        return [qqpp[2:4]; qqpp[1]]
    end

    function qp_right(q, p)
        qs, qx, qy, qz = q
        return [
            +qs -qx -qy -qz;
            +qx +qs -qz +qy;
            +qy +qz +qs -qx;
            +qz -qy +qx +qs
        ]*p
    end

    function qp_left(q, p)
        qs, qx, qy, qz = q
        return [
            +qs -qx -qy -qz;
            +qx +qs +qz -qy;
            +qy -qz +qs +qx;
            +qz +qy -qx +qs
        ]*p
    end

    qerr(a, b) = min(norm(a - b), norm(a + b))
    qp1 = q_multiplyUniversal(q, p; scalarFirst=true, right=true)
    qp2 = q_multiplyUniversal(q, p; scalarFirst=false, right=true)
    qp3 = q_multiplyUniversal(q, p; scalarFirst=true, right=false)
    qp4 = q_multiplyUniversal(q, p; scalarFirst=false, right=false)
    errv = [qerr(qp1, qp), qerr(qp2, qp), qerr(qp3, qp), qerr(qp4, qp)]

    msgs = [
        "Scalar first, right-handed (ij = k)",
        "Scalar last, right-handed (ij = k)",
        "Scalar first, left-handed (ij = -k)",
        "Scalar-last, left-handed (ij = -k)",
    ]
    println("Quaternion convention: \n"*msgs[findmin(errv)[2]])
end
