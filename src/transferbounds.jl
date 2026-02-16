@compat abstract type Space{T<:Real} end
@compat struct ChebyshevZeroOneBV{T<:Real} <: Space{T}
end
ChebyshevZeroOneBV(::Type{T}) where T<:Real = ChebyshevZeroOneBV{T}() 

@compat abstract type EntryBounds{S<:Space,T<:Real} end
@compat struct AnalyticEntryBounds{S<:Space,T<:Real}<:EntryBounds{S,T}
  sp::S
  δ::T
  P::T
  H::T
end

# Entry bounds
entrybound(a::AnalyticEntryBounds{S},j,k) where S<:ChebyshevZeroOneBV = (a.H*t(j)*exp((k-1)*a.P-(j-1)*a.δ))

# BV norm of EN
function logEnorm(a::AnalyticEntryBounds{S},NP) where S<:ChebyshevZeroOneBV
  N = NP-1; #recall we start at 0 in the polynomial basis and 1 in the matrix indexing
  log(sqrt(2)*a.H) + log(1-exp(-2a.P*N))/2 - log(1-exp(-2a.P))/2 + a.P*N-a.δ*(N+1) +
    log(sqrt(N^2 *exp(-4a.δ) + (1-2N-2N^2)*exp(-2a.δ) + (N+1)^2)/(1-exp(-2a.δ)) + 1/2) - log(1-exp(-2a.δ))/2
end

# Pointwise aliasing error
t(j) = (1+(j!=0))
aliaserror(a::AnalyticEntryBounds{S},j,k,N_interp) where S<:ChebyshevZeroOneBV =
  a.H * t(j-1) * exp(a.P*(k-1)) * 2cosh(a.δ*(j-1)) /(exp(2a.δ*N_interp)-1)



@compat abstract type LYBounds{S<:Space,T<:Real} end
@compat struct FirstOrderLYBounds{S<:Space,T<:Real}<:LYBounds{S,T}
  sp::S
  λ::T
  C1::T
end

function Snorm(b::FirstOrderLYBounds{S,T}) where {S<:ChebyshevZeroOneBV, T<:Real}
  R = 2*b.C1/(1-1/b.λ)
  ξ = exp(-R)*(1-1/b.λ)/2
  C_L = 4exp(R) * (1+R)
  C_dash =  1+b.C1/3(1-1/b.λ)
  #Γ = 1-ξ
  #Lsum_norm = (C_dash* log(C_L) + 1)/ξ
  m_lsum = ceil((2/log(b.λ)).hi);
  n_lsum = ceil((
        (4+2log(max(b.C1/(1-1/b.λ),1)*sqrt(C_L)))/ξ
      ).hi)
  Lsum_norm = (m_lsum+n_lsum)*C_dash/convert(T,3//5)
  1 + Lsum_norm*(3+C_dash)
end

# Strong convergence norm
function logSEnorm(a::EntryBounds,b::LYBounds,N)
  logenorm = logEnorm(a,N)
  logenorm - log(1/Snorm(b) - exp(logenorm))
end
