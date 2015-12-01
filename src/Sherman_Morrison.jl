function SM_solve(A, rank, row_ind::Array, b)
    n = size(A, 1)
    A0 = copy(A)
    if rank==1
        u = zeros(n,1)
        v = zeros(1,n)
        u[row_ind[1]] = 1.0
        for i in 1:n
            v[i] = A0[row_ind[1], i]
        end
        v[row_ind[1]] -= 1.0
        A0[row_ind[1], :] = 0.0
        A0[row_ind[1], row_ind[1]] = 1.0
        x1 = A0\b
        y1 = A0\u
        return x1 - (v*x1) .* y1 ./(1+v*y1)
    elseif rank==2
        u1 = zeros(n,1)
        u1[row_ind[1]] = 1.0
        v1 = zeros(1,n)
        for i in 1:n
            v1[i] = A0[row_ind[1], i]
        end
        v1[row_ind[1]] -= 1.0
        A0[row_ind[1], :] = 0.0
        A0[row_ind[1], row_ind[1]] = 1.0

        u2 = zeros(n,1)
        u2[row_ind[2]] = 1.0
        v2 = zeros(1,n)
        for i in 1:n
            v2[i] = A0[row_ind[2], i]
        end
        v2[row_ind[2]] -= 1.0
        A0[row_ind[2], :] = 0.0
        A0[row_ind[2], row_ind[2]] = 1.0

        A0_LU = lufact(A0)
        x0 = A0_LU\b
        y1 = A0_LU\u1
        z2 = A0_LU\u2
        x1 = x0 - (v1*x0) .* y1./(1+v1*y1)
        y2 = z2 - (v1*z2) .* y1./(1+v1*y1)

        return x1 - (v2*x1).*y2./(1+v2*y2)
    else
        println("rank not equal to 1 or 2 is not yet implemented!")
    end

end
