#using JuBLAS
using Base.Test

# write your own tests here
@test 1 == 1

include("../src/JuBLAS.jl")
import Base.LinAlg.gemv!

myrand{T<:Real}(::Type{T}, n...) = convert(Array{T}, rand(n...))

myrand{T<:Real}(::Type{Complex{T}}, n...) = myrand(T,n...) + im*myrand(T,n...)


tr = transpose
ctr = ctranspose

function gen_mat(trans='N', Ty=Float64, M=5, N=4)

    m = trans=='N'?M:N
    n = trans=='N'?N:M
    A = convert(Array{T,2}, randn(M,N))
    x = convert(Array{T,1}, randn(n))
    y = convert(Array{T,1}, randn(m))

    return A, x, y
end
eps1(x) = eps(x)
eps1{T<:Real}(x::Complex{T}) = eps(T)
eps1{T<:Real}(::Type{Complex{T}}) = eps(T)

# Testing Matrix vector multiply
T = Float64

α = T(1.5)
β = T(2.5)

A, x, y = gen_mat('N', T, 5, 4)
y0 = copy(y)

@test_approx_eq_eps maxabs(A*x - JuBLAS.gemv('N', 1.0, A, x)) 0.0 10*eps1(T)
@test_approx_eq_eps maxabs(α*A*x - JuBLAS.gemv('N', α, A, x)) 0.0 10*eps1(T)
@test_approx_eq_eps maxabs(α*A*x+β*y0 - JuBLAS.gemv!('N', α, A, x, β, y)) 0.0 10*eps1(T)

A, x, y = gen_mat('T', T, 5, 4)
y0 = copy(y)


@test_approx_eq_eps maxabs(tr(A)*x - JuBLAS.gemv('T', 1.0, A, x)) 0.0 10*eps1(T)
@test_approx_eq_eps maxabs(α*tr(A)*x - JuBLAS.gemv('T', α, A, x)) 0.0 10*eps1(T)
@test_approx_eq_eps maxabs(α*tr(A)*x+β*y0 - JuBLAS.gemv!('T', α, A, x, β, y)) 0.0 10*eps1(T)




T = Complex{Float64}

α = T(1.5)
β = T(2.5)

A, x, y = gen_mat('N', T, 5, 4)
y0 = copy(y)


@test_approx_eq_eps maxabs(A*x - JuBLAS.gemv('N', 1.0, A, x)) 0.0 10*eps1(T)
@test_approx_eq_eps maxabs(α*A*x - JuBLAS.gemv('N', α, A, x)) 0.0 10*eps1(T)
@test_approx_eq_eps maxabs(α*A*x+β*y0 - JuBLAS.gemv!('N', α, A, x, β, y)) 0.0 10*eps1(T)


A, x, y = gen_mat('T', T, 5, 4)
y0 = copy(y)

@test_approx_eq_eps maxabs(tr(A)*x - JuBLAS.gemv('T', 1.0, A, x)) 0.0 10*eps1(T)
@test_approx_eq_eps maxabs(α*tr(A)*x - JuBLAS.gemv('T', α, A, x)) 0.0 10*eps1(T)
@test_approx_eq_eps maxabs(α*tr(A)*x+β*y0 - JuBLAS.gemv!('T', α, A, x, β, y)) 0.0 10*eps1(T)


A, x, y = gen_mat('T', T, 5, 4)
y0 = copy(y)

@test_approx_eq_eps maxabs(ctr(A)*x - JuBLAS.gemv('C', 1.0, A, x)) 0.0 10*eps1(T)
@test_approx_eq_eps maxabs(α*ctr(A)*x - JuBLAS.gemv('C', α, A, x)) 0.0 10*eps1(T)
@test_approx_eq_eps maxabs(α*ctr(A)*x+β*y0 - JuBLAS.gemv!('C', α, A, x, β, y)) 0.0 10*eps1(T)

#=
=#
