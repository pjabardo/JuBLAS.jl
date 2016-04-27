

function gemm!{T<:Number}(transA::Char, transB::Char, alpha, A::AbstractMatrix{T}, 
                          B::AbstractMatrix{T}, beta, C::AbstractMatrix{T})
    transA = uppercase(transA)
    transB = uppercase(transB)

    M = size(A,1)
    N = size(B,2)
    K = size(A,2)
    
    notA = transA == 'N'
    notB = transB == 'N'
    conjA = transA == 'C'
    conjB = transB == 'C'
    
    
    if notA
        nrowA = M
        ncolA = K
    else
        nrowA = K
        ncolA = M
    end

    if notB
        nrowB = K
    else
        nrowB = N
    end
    
    if M==0 || N==0 || ((alpha==0 || K==0) && beta==1)
        return C
    end

    if alpha==0
        if beta==0
            for j = 1:N
                for i = 1:M
                    C[i,j] = zero(T)
                end
            end
        else
            for j = 1:N
                for i = 1:M
                    C[i,j] = beta*C[i,j]
                end
            end
        end
        return C
    end

    if notB
        if notA
            for j = 1:N
                if beta==0
                    for i = 1:M
                        C[i,j] = zero(T)
                    end
                elseif beta != 1
                    for i = 1:M
                        C[i,j] = beta*C[i,j]
                    end
                end
                for l = 1:K
                    if B[l,j] != zero(T)
                        temp = alpha * B[l,j]
                        for i = 1:M
                            C[i,j] += temp*A[i,l]
                        end
                    end
                end
            end
        elseif conjA
            for j = 1:N
                for i = 1:M
                    temp = zero(T)
                    for l = 1:K
                        temp += conj(A[l,i]) * B[l,j]
                    end
                    if beta==0
                        C[i,j] = alpha*temp
                    else
                        C[i,j] = alpha*temp + beta*C[i,j]
                    end
                end
            end
        else
            for j = 1:N
                for i = 1:M
                    temp = zero(T)
                    for l = 1:K
                        temp += A[l,i] * B[l,j]
                    end
                    if beta==0
                        C[i,j] = alpha*temp
                    else
                        C[i,j] = alpha*temp + beta*C[i,j]
                    end
                end
            end
        end
    elseif notA
        if conjB
            for j = 1:N
                if beta==0
                    for i = 1:M
                        C[i,j] = zero(T)
                    end
                elseif beta!=1
                    for i = 1:M
                        C[i,j] = beta*C[i,j]
                    end
                end
                for l = 1:K
                    if B[j,l] != zero(T)
                        temp = alpha * conj(B[j,l])
                        for i = 1:M
                            C[i,j] += temp*A[i,l]
                        end
                    end
                end
            end
        else
            for j = 1:N
                if beta==0
                    for i = 1:M
                        C[i,j] = zero(T)
                    end
                elseif beta != 1
                    for i = 1:M
                        C[i,j] = beta*C[i,j]
                    end
                end
                for l = 1:K
                    if B[j,l] != zero(T)
                        temp = alpha*B[j,l]
                        for i = 1:M
                            C[i,j] += temp*A[i,l]
                        end
                    end
                end
            end
        end
    elseif conjA
        if conjB
            for j = 1:N
                for i = 1:M
                    temp = zero(T)
                    for l = 1:K
                        temp += conj(A[l,i]) * conj(B[j,l])
                    end
                    if beta == 0
                        C[i,j] = alpha*temp
                    else
                        C[i,j] = alpha*temp + beta*C[i,j]
                    end
                end
            end
        else
            for j = 1:N
                for i = 1:M
                    temp = zero(T)
                    for l = 1:K
                        temp += conj(A[l,i])*B[j,l]
                    end
                    if beta == 0
                        C[i,j] = alpha*temp
                    else
                        C[i,j] = alpha*temp + beta*C[i,j]
                    end
                end
            end
        end
    else
        if conjB
            for j = 1:N
                for i = 1:M
                    temp = zero(T)
                    for l = 1:K
                        temp += A[l,j]*conj(B[j,l])
                    end
                    if beta==0
                        C[i,j] = alpha*temp
                    else
                        C[i,j] = alpha*temp + beta*C[i,j]
                    end
                end
            end
        else
            for j = 1:N
                for i = 1:M
                    temp = zero(T)
                    for l = 1:K
                        temp += A[l,i]*B[j,l]
                    end
                    if beta==0
                        C[i,j] = alpha*temp
                    else
                        C[i,j] = alpha*temp + beta*C[i,j]
                    end
                end
            end
        end
    end
    return C


end

function gemm{T<:Number}(transA::Char, transB::Char, alpha, A::AbstractMatrix{T}, 
                         B::AbstractMatrix{T})
    gemm!(transA, transB, alpha, A, B, zero(T),
          similar(A, T, size(A, transA=='N'?1:2),
                  size(B, transB=='N'?2:1)))
end

