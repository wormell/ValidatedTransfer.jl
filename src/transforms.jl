# from FastTransforms
function backend_fft_pow2!{T}(x::Vector{T})
    n,big2,bigpi=length(x),2one(T),T(pi)
    nn,j=n÷2,1
    @inbounds for i=1:2:n-1
        if j>i
            x[j], x[i] = x[i], x[j]
            x[j+1], x[i+1] = x[i+1], x[j+1]
        end
        m = nn
        while m ≥ 2 && j > m
            j -= m
            m = m÷2
        end
        j += m
    end
    logn = 2
    while logn < n
        piθ=-bigpi*big2/logn
        wtemp = sin(piθ/2)
        wpr, wpi = -2wtemp^2, sin(piθ)
        wr, wi = one(T), zero(T)
        for m=1:2:logn-1
            @inbounds for i=m:2logn:n
                j=i+logn
                mixr, mixi = wr*x[j]-wi*x[j+1], wr*x[j+1]+wi*x[j]
                x[j], x[j+1] = x[i]-mixr, x[i+1]-mixi
                x[i], x[i+1] = x[i]+mixr, x[i+1]+mixi
            end
            wr = (wtemp=wr)*wpr-wi*wpi+wr
            wi = wi*wpr+wtemp*wpi+wi
        end
        logn = logn << 1
    end
    return x
end

function fft_pow2{T}(x::Vector{T})
    y = FastTransforms.interlace(real(x),imag(x))
    backend_fft_pow2!(y)
@compat    return complex.(y[1:2:end],y[2:2:end])
end

## TAYLOR

taylorinterpolationerror(r,Ninterp::Int,M,R) = M*(r/R)^Ninterp/(1-r/R)

function taylortransform{T}(vals::Vector{Complex{T}},r::T,M::T,R::T)
  Ninterp = length(vals)
  cfs = fft_pow2(vals)/Ninterp
  err= erroradjustment(Complex{T},taylorinterpolationerror(r,Ninterp,M,R))
  powr = one(T)
  @inbounds for i = 1:Ninterp
    cfs[i] *= powr
    cfs[i] += err
    err /= R
    powr /= r
  end
  cfs
end

## CHEBYSHEV

function chebyshevtransform{T}(a::AbstractArray{Complex{T}})
  N = length(a) #big(length(a))
    ispow2(N) || error("N not a power of 2")
    c = fft_pow2([a; flipdim(a,1)])
    d = c[1:N]
 @compat    d .*= exp.((-im*convert(T,pi)).*(0:N-1)./(2*N))
     d[1] = d[1] / 2 #sqrt(big(2))
     scale!(inv(convert(T,N)), d)
end
@compat chebyshevtransform{T<:Real}(a::AbstractArray{T}) = real(chebyshevtransform(complex(a)))


function cos_vec!(A::Vector)
  for i in eachindex(A)
    A[i] = cos(A[i])
  end
  A
end
chebyshevpoints(N_interp,T) = cos_vec!(collect(1:2:2N_interp-1)*(T(pi)/2N_interp));

chebyshev_shifted_sq(x,k) = (-1)^k * cos(2k*asin(x))
chebyshev_shifted(x,k) = chebyshev_shifted_sq(sqrt(x),k)

function chebyshev_shifted_by_squaring{T}(x::T,k)
  tot=k;
  s0sq = 4x*(1-x)
  cn = 2x-1; sn_div = one(T)
  ctot = one(T); stot_div = zero(T)
  while tot>0
    if isodd(tot)
      ctot,stot_div = ctot*cn - stot_div*sn_div*s0sq, stot_div*cn - ctot*sn_div
    end
    sn_div *= 2cn; cn *= 2cn; cn-=1;
    tot ÷=2
  end
  ctot
end

function chebyshev_shifted_multiel{T}(x::T,K)
  Tk = Array{T}(K)
  mul = 4x-2
  K > 0 && (Tk[1] = one(typeof(x)))
  K > 1 && (Tk[2] = 2x-1)
  @inbounds for kk = 1:div(K-1,2)
    Tk[2kk+1] = 2Tk[kk+1]^2-1
    Tk[2kk+2] = mul*Tk[2kk+1]-Tk[2kk]
  end
  isodd(K) && (Tk[K] = 2Tk[div(K,2)+1]^2-1)
  Tk
end

function fejerweights(n,T::DataType=Float64)
  v = Array{Complex{T}}(n)
  v[1] = 2
  # twid = exp(-im*T(pi)/n) # conjugating here to get ifft
  # twidp = copy(twid)
  @inbounds for k = 1:div(n-1,2)
    twidp = exp(-im*k*T(pi)/n)
    v[k+1] = 2twidp/(1-4k^2)
    # twidp *= twid
  end
  mod(n,2)==0 && (v[div(n,2)+1]=0)
  @inbounds for k = 1:div(n-1,2)
    v[n-k+1] = conj(v[k+1])
  end
   real(fft_pow2(v))/n
end
 # fejerweights{T<:Real}(n,::Type{Complex{T}}) = fejerweights(n,T)

fejererror(N,M,ρ) = 2M/(ρ-1)*ρ^(-N)
