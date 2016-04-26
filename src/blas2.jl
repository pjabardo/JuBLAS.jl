


function gemv!{T<:Number}(trans::Char, alpha, A::AbstractMatrix{T},
                          x::AbstractVector{T}, beta, y::AbstractVector{T})

    trans = uppercase(trans)
    o = one(T)
    z = zero(T)
    info = 0
    M = size(A,1)
    N = size(A,2)
    
    if !(trans=='N') && !(trans=='T') && !(trans=='C')
        info = 1
    elseif M < 0
        info = 2
    elseif N < 0
        info = 3
    end

    noconj = false
    if trans=='T'
        noconj = true
    end
    
    
    if M==0 || N==0 || (alpha==z && beta==o)
        return y
    end

    if trans=='N'
        lenx = N
        leny = M
    else
        lenx = M
        leny = N
    end

    kx = 1
    ky = 1
    
    if beta != 1
        if beta == 0
            for i in 1:leny
                y[i] = zero(T)
            end
        else
            for i = 1:leny
                y[i] = beta*y[i]
            end
        end
    end
    if alpha==z
        return y
    end
    
    if trans=='N'

        for j = 1:N
            if x[j] != z
                temp = alpha * x[j]
                for i = 1:M
                    y[i] +=  temp*A[i,j]
                end
            end
        end

    else

        for j = 1:N
            temp = z
            if noconj
                for i = 1:M
                    temp += A[i,j]*x[i]
                end
            else
                for i = 1:M
                    temp += conj(A[i,j]) * x[i]
                end
            end
            y[j] += alpha*temp
        end
    end

    return y
    
end

function gemv{T<:Number}(trans::Char, alpha, A::AbstractMatrix{T},x::AbstractVector{T})
    gemv!(trans, one(T), A, x, zero(T), similar(x, T, size(A, (trans=='N'?1:2))))
end


