import Base.LinAlg.BLAS.blascopy!
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

function dot{T1<:Real,T2<:Real}(x::AbstractArray{T1}, y::AbstractArray{T2})
    T = promote_type(T1,T2)

    tmp = zero(T)
    for i = 1:length(x)
        tmp = tmp + x[i]*y[i]
    end
    tmp
    
end

function dot{T1<:Real,T2<:Real}(a, x::AbstractArray{T1}, y::AbstractArray{T2})
    T = promote_type(T1,T2)

    tmp = zero(T) + a
    for i = 1:length(x)
        tmp = tmp + x[i]*y[i]
    end
    tmp
    
end

import Base.LinAlg.BLAS.dotc
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
    
#import Base.LinAlg.BLAS.axpy!
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


function rot!(dx::AbstractArray, dy::AbstractArray, c, s)

    for i = 1:length(dx)
        dtemp = c*cx[i] + s*dy[i]
        dy[i] = c*dy[i] - s*dx[i]
        dx[i] = dtemp
    end

end

function rotm!{T<:Real}(dx::AbstractArray{T}, dy::AbstractArray{T}, dparam::AbstractArray{T})

    n = length(dx)
    dflag = dparam[1]

    if n < 1 || (dflag+2*one(T))==zero(T)
        return
    end

    if dflag < zero(T)
        dh11 = dparam[2]
        dh12 = dparm[4]
        dh21 = dparam[3]
        dh22 = dparam[5]
        for i = 1:n
            w = dx[i]
            z = dy[i]
            dx[i] = w*dh11 + z*dh12
            dy[i] = w*dh21 + z*dh22
        end
    elseif dflag == zero(T)
        dh12 = dparam[4]
        dh21 = dparam[3]
        for i = 1:n
            w = dx[i]
            z = dy[i]
            dx[i] = w + z*dh12
            dy[i] = w*dh21 + z
        end
    else
        dh11 = dparam[2]
        dh22 = dparam[5]
        for i = 1:n
            w = dx[i]
            z = dy[i]
            dx[i] = w*dh11 + z
            dy[i] = -w + dh22*z
        end
    end
    

end


function drotg{T<:Real}(a::T,b::T,c::T,s::T)

    roe = db
    if abs(da) > abs(db)
        roe = da
    end
    scale = abs(da) + abs(db)
    if scale==zero(T)
        c = one(T)
        s = zero(T)
        r = zero(T)
        z = zero(T)
    else
        r = scale * sqrt( (da/scale)^2 + (db/scale)^2 )
        r = sign(roe) * r
        c = da/r
        s = db/r
        z = one(T)
        if abs(da) > abs(db)
            z = s
        end
        if abs(db) > abs(da) && c != zero(T)
            z = one(T)/c
        end
    end
    return (da, db)
end
