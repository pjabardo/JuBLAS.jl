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

function dot(x::AbstractArray{T1}, y::AbstractArray{T2}) where {T2<:Real,T1<:Real}
    T = promote_type(T1,T2)

    tmp = zero(T)
    for i = 1:length(x)
        tmp = tmp + x[i]*y[i]
    end
    tmp
    
end

function dot(a, x::AbstractArray{T1}, y::AbstractArray{T2}) where {T2<:Real,T1<:Real}
    T = promote_type(T1,T2)

    tmp = zero(T) + a
    for i = 1:length(x)
        tmp = tmp + x[i]*y[i]
    end
    tmp
    
end

import Base.LinAlg.BLAS.dotc
function dotc(x::AbstractArray{Complex{T}}, y::AbstractArray{Complex{T}}) where {T<:Real}
    tmp = zero(Complex{T})
    for i = 1:length(x)
        tmp = tmp + conj(x[i]) * y[i]
    end
    tmp
end


import Base.LinAlg.BLAS.dotu
function dotu(x::AbstractArray{Complex{T}}, y::AbstractArray{Complex{T}}) where {T<:Real}
    tmp = zero(Complex{T})
    for i = 1:length(x)
        tmp = tmp + x[i] * y[i]
    end
    tmp
end

    
import Base.LinAlg.BLAS.nrm2
function nrm2(x::AbstractArray{T}) where {T<:Number}
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
function asum(x::AbstractArray{T}) where {T<:Number}

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

function rotm!(dx::AbstractArray{T}, dy::AbstractArray{T}, dparam::AbstractArray{T}) where {T<:Real}

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


function rotg(da::T,db::T,c::T,s::T) where {T<:Real}

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



function rotmg!(dd1::T, dd2::T, dx1::T, dy1::T, dparam::AbstractVector{T}) where {T<:Real}

    two = 2*one(T)
    gam = convert(T, 4096)
    gamsq = convert(T, 16777216)
    rgamsq = convert(T, 5.9604645e-8)

    if dd1 < z
        dflag = -one(T)
        dh11 = zero(T)
        dh12 = zero(T)
        dh21 = zero(T)
        dh22 = zero(T)

        dd1 = zero(T)
        dd2 = zero(T)
        dx1 = zero(T)
    else
        dp2 = dd2 * dy1
        if dp2==zero(T)
            dflag = -two
            dparam[1] = dflag
            return dd1, dd2, dx1, dy1, dparam
        end
        dp1 = dd1*dx1
        dq2 = dp2*dy1
        dq1 = dp1*dx1
        
        if abs(dq1) > abs(dq2)
            dh21 = -dy1/dx1
            dh12 = dp2 / dp1
            
            du = one(T) - dh12*dh21
            if du > zero(T)
                dflag = zero(T)
                dd1 = dd1/du
                dd2 = dd2/du
                dx1 = dx1*du
            end
        else
            if dq2 < zero(T)
                dflag = -one(T)
                dh11 = zero(T)
                dh11 = zero(T)
                dh12 = zero(T)
                dh21 = zero(T)
                dh22 = zero(T)
                
                dd1 = zero(T)
                dd2 = zero(T)
                dx1 = zero(T)
            else
                dflag = one(T)
                dh11 = dp1/dp2
                dh22 = dx1 / dy1
                du = one(T) + dh11*dh22
                dtemp = dd2/du
                dd2 = dd1/du
                dd1 = dtemp
                dx1 = dy1*du
            end
        end
        
        if dd1 != zero(T)
            while (dd1 <= rgamsq || dd1 >= gamsq)
                if dflag==zero(T)
                    dh11 = one(T)
                    dh22 = one(T)
                    dflag = -one(T)
                else
                    dh21 = -one(T)
                    dh12 = one(T)
                    dflag = -one(T)
                end
                if dd1 <= rgamsq
                    dd1 = dd1 * gam*gam
                    dx1 = dx1/gam
                    dh11 = dh11/gam
                    dh12 = dh12/gam
                else
                    dd1 = dd1 / (gam*gam)
                    dx1 = dx1 * gam
                    dh11 = dh11*gam
                    dh12 = dh12*gam
                end
                
            end
        end

        if dd2 != zero(T)
            while (abs(dd2) <= rgamsq || abs(dd2) >= gamsq)
                if dflag==zero(T)
                    dh11 = one(T)
                    dh22 = one(T)
                    dflag = -one(T)
                else
                    dh21 = -one(T)
                    dh12 = one(T)
                    dflag = -one(T)
                end
                if abs(dd2) <= rgamsq
                    dd2 = dd2*gam*gam
                    dh21 = dh21 / gam
                    dh22 = dh22 / gam
                else
                    dd2 = dd2 / (gam*gam)
                    dh21 = dh21*gam
                    dh22 = dh22*gam
                end
            end
        end

    end

    if dflag < zero(T)
        dparam[2] = dh11
        dparam[3] = dh21
        dparam[4] = dh12
        dparam[5] = dh22
    elseif dflag==zero(T)
        dparam[3] = dh21
        dparam[4] = dh12
    else
        dparam[2] = dh11
        dparam[5] = dh22
    end

    dparam[1] = dflag

    return dd1, dd2, dx1, dy1, dparam
               
    
end



function swap!(x::AbstractArray, y::AbstractArray)

    n = length(x)
    for i = 1:n
        dtemp = x[i]
        dx[i] = dy[i]
        dy[i] = dtemp
    end
    
end

