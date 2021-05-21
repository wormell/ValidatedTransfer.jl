 addprocs(3)
@everywhere using Base.Test, Compat, IntervalArithmetic, ValidatedTransfer

# module Gauss
#   using IntervalArithmetic, ValidatedTransfer, Compat
#   # println(names(IntervalArithmetic))
#   # println(names(ValidatedTransfer))
#
#   λ = @biginterval 2
#   λ̌ = sqrt(λ)
#   C1 = @biginterval log(2)
#   δ = @biginterval acosh(16/5)
#   P = @biginterval acosh(2log2(1+2^1.1)-1)
#   H = @biginterval 2
#
#   sp = ChebyshevZeroOneBV(BigInterval)
#   entrybounds = AnalyticEntryBounds(sp,δ,P,H)
#   lybounds = FirstOrderLYBounds(sp,λ,C1)
#
#
#   # const biglog2 = (@biginterval(log(2)))::BigInterval
#   Compat.@compat immutable g end; Compat.@compat (::g){T}(x::T) = (1/(exp2(x)-1))::T
#   Compat.@compat immutable dg end; Compat.@compat (::dg){T}(x::T) = begin e2x = exp2(x); e2x1 = e2x-1; (-log(T(2))*e2x/(e2x1*e2x1))::T end
#   Compat.@compat immutable γa end; Compat.@compat (::γa){T}(x::T) = (log2(1+1/x))::T
#   Compat.@compat immutable h end; Compat.@compat (::h){T}(x::T) = (exp2(x)-1)::T
#   Compat.@compat immutable dh end; Compat.@compat (::dh){T}(x::T) = (log(T(2))*exp2(x))::T
#   α = IntervalArithmetic.@biginterval 1;
#   α1 = IntervalArithmetic.@biginterval -1;
#   G0 = IntervalArithmetic.@biginterval 1/log(2);
#   G1 = IntervalArithmetic.@biginterval 0.1
#   saferad = IntervalArithmetic.@biginterval 0.4;


# LANFORD TEST

@everywhere @eval begin
  N = 100
  N_interp = nextpow2(N)

  prec = 8round(Int,(N+512)*0.2/8)
  @compat setprecision(BigFloat,prec)

  import Base.sinpi
  Base.sinpi(x::IntervalArithmetic.Interval) = sin(pi*x)

  ## MAP CONSTANTS
  @compat lan_v(x,b) = (5 - sqrt.(25-8(x+b)))/2
  @compat lan_dv(x,b) = 2./sqrt.(25-8(x+b))
  λ = @biginterval 3/2
  λ̌ = sqrt(λ)
  C1 = @biginterval 4/9
  δ = @biginterval acosh( 7/4)
  P = @biginterval acosh(4 - sqrt(6))
  H = @biginterval sqrt((7+sqrt(33))/2) #5/2sqrt(39-24sqrt(2)) + 3/2sqrt(40sqrt(2)-55)

  # Pointwise aliasing error
  t(j) = (1+(j!=0))
  aliaserror(j,k,N) = (H * t(j-1) * exp(P*(k-1)) * 2cosh(δ*(j-1)) /(exp(2δ*N)-1)).hi
  entrybound(j,k) = (H*t(j)*exp((k-1)*P-(j-1)*δ)).hi

  @compat @assert precision(BigFloat) == prec

  function acos_lan_v_norm{T}(A::Vector{T},b)
    B = similar(A)
    for i in eachindex(B)
      B[i] = acos(2lan_v((A[i]+1)/2,b)-1)
    end
    B
  end

  function lan_dv_norm{T}(A::Vector{T},b)
    B = similar(A)
    for i in eachindex(B)
      B[i] = lan_dv((A[i]+1)/2,b)
    end
    B
  end

  cheb_pts = ValidatedTransfer.chebyshevpoints(N_interp,BigInterval)
  acos_v0_cheb_pts = acos_lan_v_norm(cheb_pts,0); acos_v1_cheb_pts = acos_lan_v_norm(cheb_pts,1);
  dv0_cheb_pts = lan_dv_norm(cheb_pts,0); dv1_cheb_pts = lan_dv_norm(cheb_pts,1);
  helperdata = [acos_v0_cheb_pts,dv0_cheb_pts,acos_v1_cheb_pts,dv1_cheb_pts]

  println(size(acos_v0_cheb_pts))
  lanford_transferfn(i,k,helperdata) = cos((k-1)*helperdata[1][i])*helperdata[2][i] + cos((k-1)*helperdata[3][i])*helperdata[4][i]
end

@time M = maketransfer(N,lanford_transferfn,helperdata,aliaserror,entrybound,workers(),N_interp)
@time S = getsolution(M,workers())
@time savezip(tempname(),S)
