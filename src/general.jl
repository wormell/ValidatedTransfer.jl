const BigInterval = IntervalArithmetic.Interval{BigFloat} #typealias BigInterval IntervalArithmetic.Interval{BigFloat}

function upbound(x)
  xu = _upbound(x)
  assert(xu≥0)
  xu
end
_upbound(x::Real) = x
_upbound(x::Interval) = x.hi

roundupbound(x) = ceil(Int,upbound(x))

prec{T}(::Type{Interval{T}}) = eps(T)
prec(T) = eps(T)

erroradjustment(T,err) = zero(T)
erroradjustment{T}(::Type{Interval{T}},err) = Interval(-upbound(err),upbound(err))
erroradjustment{IT<:Interval}(::Type{Complex{IT}},err) =
  complex(erroradjustment(IT,err),erroradjustment(IT,err))
erroradjustment(err) = erroradjustment(typeof(err),err)

intervalsign(S) =  (S > 0) ? 1 : (S < 0) ? -1 : error("$S contains zero")

# move to IntervalArithmetic
# import Base.log, Base.hypot
import Base.exp2, Base.sqrt

# VERSION < v"0.5" && (Base.hypot(x::Interval,y::Interval) = sqrt(x*x+y*y))
function Base.sqrt{T<:Interval}(z::Complex{T})
        (x,y) = reim(z)
        nz = hypot(x,y)
        if 0 ∈ y && x.lo < 0
          error("Not implemented")
          # if y == 0
          #   return 0∈x complex(y,sqrt(-x))
          # (xlh,ylh) = reim(csqrt(Interval(complex(x.lo,y.hi))))
          # (xhh,yhh) = reim(csqrt(Interval(complex(x.hi,y.hi))))
          # if y.lo == 0
          #   re = Interval(0,xhh.hi)
          #   im = Interval(sqrt(Interval(-x.hi)),xlh.hi)
          # else
          #   (xll,yll) = reim(csqrt(Interval(complex(x.lo,y.lo))))
          #   (xhl,yhl) = reim(csqrt(Interval(complex(x.hi,y.lo))))
          #   re = Interval(0,max(xhl.hi,xhh.hi))
          #   im = Interval(yll.lo,ylh.hi)
          # end
        else
          re = sqrt((nz+x)/2)
          im = y/2re
        end
        complex(re,im)
end
#
# Base.log{T<:Interval}(z::Complex{T}) = complex(log(abs(z)),atan2(imag(z),real(z)))
#
function Base.exp2{T<:Real}(z::Complex{T})
    er = exp2(real(z))
    theta = imag(z) * log(convert(T, 2))
    Complex(er*cos(theta), er*sin(theta))
end


# Akiyama–Tanigawa algorithm for second Bernoulli numbers B+n
function bernoulli2k(p::Integer,T)
  B2 = Array{T}(p+1)
  A = Array{BigInt}(2p+1)
  fact = one(BigInt)
  for m = 0:2p
      A[m+1] = fact
      fact *= m+1
      for j = m:(-1):1
          A[j] = j*((m+1)*A[j]-A[j+1])
      end
      iseven(m) && (B2[div(m,2)+1] = T(A[1])/fact
      )
  end
  return B2
end
