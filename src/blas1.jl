import Base.LinAlg.BLAS.blascopy
function blascopy!(x::AbstractArray, y::AbstractArray)
    for i = 1:length(x)
        y[i] = x[i]
    end
    dy
end

import Base.LinAlg.BLAS.scal!
function scal!(a, x::AbstractArray)
    for i = 1:length(x)
        x[i] = a * x[i]
    end
end

function dot{T<:Real}(x::AbstractArray{T}, y::AbstractArray{T})


    tmp = zero(T)
    for i = 1:length(x)
        tmp = tmp + x[i]*y[i]
    end
    tmp
    
end

import Base.LinAlg.BLAS.doct
function dotc{T<:Real}(x::AbstractArray{Complex{T}}, y::AbstractArray{Complex{T}})
    tmp = zero(Complex{T})
    for i = 1:length(x)
        tmp = tmp + conj(x[i]) * y[i]
    end
    tmp
end


import Base.LinAlg.BLAS.dotu
function dotu{T<:Real}(x::AbstractArray{Complex{T}}, y::AbstractArray{Complex{T}})
    tmp = zero(Complex{T})
    for i = 1:length(x)
        tmp = tmp + x[i] * y[i]
    end
    tmp
end

    
import Base.LinAlg.BLAS.nrm2
function nrm2{T<:Number}(x::AbstractArray{T})
    n = length(x)
    if n < 1
        return abs(zero(T))
    elseif n == 1
        return abs(x[1])
    else
        tmp = abs(x[1])
        for i = 2:n
            tmp = tmp + abs(x[i])^2
        end
        return sqrt(tmp)
    end
end

    
import Base.LinAlg.BLAS.asum
function asum{T<:Number}(x::AbstractArray{T})

    n = length(x)
    if n < 1
        return abs(zero(T))
    elseif n == 1
        return abs(x[1])
    else
        tmp = abs(x[1])
        for i = 2:n
            tmp = tmp + abs(x[i])
        end
        return tmp
    end
    
end
    
import Base.LinAlg.BLAS.axpy!
function axpy!(a, x::AbstractArray, y::AbstractArray)

    for i = 1:length(x)
        y[i] = y[i] + a*x[i]
    end
    dy
        
end

import Base.LinAlg.BLAS.iamax
function iamax(x::AbstractArray)
    n = length(x)
    if n < 1
        return 0
    elseif n == 1
        return 1
    end

    idx = 1
    xmax = abs(x[1])
    for i = 2:n
        xai = abs(x[i])
        if xai > xmax
            idx = i
            xmax = xmax
        end
    end
    idx
end


