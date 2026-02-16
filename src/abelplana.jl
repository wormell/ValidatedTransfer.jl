## ALL VERY DEPRECATED

struct AbelPlanaPlan{T,Tg,Tdg,Tγ,Th,Tdh}
  g::Tg
  dg::Tdg
  γ::Tγ
  h::Th
  dh::Tdh
  α::T
  α1::T
  G0::T
  G1::T
  nmin::Int
  saferad::T
  Θmax::T
  sbdmin::T
end

AbelPlanaPlan(args...;kwargs...) = _AbelPlanaPlan(args...;kwargs...)
function _AbelPlanaPlan(g,dg,γ,h,dh,α,α1,G0,G1;nmin::Int=0,saferad=1,Θmax=2)
  (α,α1,G0,G1,saferad,Θmax) = promote(α,α1,G0,G1,saferad,Θmax)
  T = typeof(α)
  sbdmin = max(-log(saferad),-log(γ(T(nmin))),log(T(4)))
  G1 != 0 && (sbdmin = max(sbdmin,log((3T(pi)/4sin(3T(pi)/8)*G1*abs(α1))/(G0*α))/(α-α1)))
  Θmax = min(Θmax,T(pi)/2α,T(pi)/3)
  AbelPlanaPlan{T,typeof(g),typeof(dg),typeof(γ),typeof(h),typeof(dh)}(g,dg,γ,h,dh,α,α1,G0,G1,convert(Int,nmin),saferad,Θmax,sbdmin)
end

function transfer_basisfun_parameters{T}(ap::AbelPlanaPlan{T},x::T,k::Integer,tol=200prec(T),debug::Bool=false)
  α=ap.α; α1 = ap.α1; G0 = ap.G0; G1 = ap.G1

  hx = ap.h(x); bigpi = T(pi)


  L = max(one(T),T(k)/40) #one(T) #whatever
  sbd = ap.sbdmin#max(ap.sbdmin,-log(L*ap.α*ap.G0/1024k)/(α+one(T)/2))# lol that constant
  exp_msbd_sqrt = exp(-sbd/2)
  Θ = min(asin(L/max(k,1)*(cos(ap.Θmax)/exp_msbd_sqrt-exp_msbd_sqrt)),
        ap.Θmax)
  Δ = Θ/3; θ = Θ-Δ

  swidth = min(sqrt((2Δ/θ+1)^2 - 1),1//4)
  # TODO: choose sminus szero splus so that the width of R1 is
  # comparable with length re: Bernstein ellipse
  mminus = roundupbound(ap.g(exp(-sbd))); sminus = -log(ap.γ(mminus+one(T)/2+hx))
  mzero = roundupbound(ap.g(exp(-swidth-sminus))-hx); szero = -log(ap.γ(mzero+one(T)/2+hx))
  mplus = roundupbound(ap.g(exp(-swidth-szero))-hx); splus = -log(ap.γ(mplus+one(T)/2+hx))

  debug && println((mminus,mzero,mplus,sminus,szero,splus))
  debug && println("sbd=",sbd," big Θ=",Θ)
  ϕmax = exp(L)
  G0pert = G1*α/abs(α1) * exp(-(α-α1)*sbd)
  dtilde = (G0 - α*(θ+2Δ)/sin(α*(θ+2Δ)/2)*G0pert)*
    min(exp(α*sminus)*tan(α*Δ),1/2(G0+G0pert))

  debug && println("ϕmax=",ϕmax," G0pert=",G0pert," dtilde=",dtilde)
  # C1 INTEGRAL
  semiminor1 = min((splus-szero)/θ,(szero-sminus)/θ)
  ρ1 = min(2+sqrt(T(5))/4,semiminor1+sqrt(semiminor1^2+1))
  M1 = ϕmax*( exp(-sbd)*(coth(π*dtilde)+1)/2 +
        (mplus-mminus)/(2bigpi*dtilde) * exp(-splus) *
          (α*G0*exp(α*splus)+abs(α1)*G1*exp(α1*splus)) /
          (α*G0*exp(α*sbd)-abs(α1)*G1*exp(α1*sbd))
      )
      # println(M1,ρ1)
      # println("Unsafe?: ",-log(4M1*tol/3/(ρ1-1))/log(ρ1))
  N1 = nextpow2(roundupbound(-log(tol/3*(ρ1-1)/2M1)/log(ρ1)))
  debug && println("ρ1=",ρ1," M1=",M1," N1=",N1)

  # C2 INTEGRAL
  Kθ = sin(α*θ)*(G0-G0pert)
  xstar = max(log(-log(tol/3*Kθ*α/ϕmax)/Kθ)/α,szero)
  #  exp_mxst = exp(-xstar-θ*im)
  #  -exp_mxst.*chebyshev_shifted.(exp_mxst,k)./expm1.(-2bigpi*im*(ap.g.(exp_mxst) - hx))|>println
  Dx = 10Δ; nDx = roundupbound((xstar-szero)/2Dx)
  debug && println("xstar=",xstar," Dx=",Dx," nDx=",nDx)
  # if nDx > 0
    semimajor2 = 1 + (szero-sminus)/Dx; semiminor2 = Δ/Dx
    ρ2 = min(semimajor2+sqrt(semimajor2^2-1),semiminor2+sqrt(semiminor2^2+1))
    M2 = ϕmax*exp(-sminus)/(exp(sin(α*(θ-Δ))*(G0-G0pert)*exp(α*sminus))-1)
    N2 = nextpow2(roundupbound(-log(tol/3nDx*(ρ2-1)/2M2)/log(ρ2)))
    debug && println("ρ2=",ρ2," M2=",M2," N2=",N2)
  # end

  exp_msend = exp(-szero-(2nDx+1)*Dx-im*θ)
  I2a2err = abs(exp_msend)*ϕmax/abs(expm1(-2bigpi*im*(ap.g(exp_msend)-hx)))
  debug && println("I2a2err=",I2a2err)

  α,α1,G0,G1,hx,bigpi,θ,mminus,mzero,mplus,szero,
    ρ1,M1,N1,ρ2,M2,N2,Dx,nDx,I2a2err
end


function transfer_basisfun{T}(ap::AbelPlanaPlan{T},x::T,k::Integer,tol=200prec(T);debug::Bool=false)
  @assert x>=0 k>=0

  α,α1,G0,G1,hx,bigpi,θ,mminus,mzero,mplus,szero,
    ρ1,M1,N1,ρ2,M2,N2,Dx,nDx,I2a2err = transfer_basisfun_parameters(ap,x,k,tol,debug)

  # C1 INTEGRAL
  chebpts = chebyshevpoints(N1,T)
  exp_ms = Array{Complex{T}}(undef,N1)
  I2a1_vals = Array{Complex{T}}(undef,N1)
  g_exp_ms_h = Array{Complex{T}}(undef,N1)
  for i = eachindex(exp_ms)
    exp_ms[i] = exp(-szero - im*θ*(1+chebpts[i])/2)
    g_exp_ms_h[i] = ap.g(exp_ms[i]) - hx
    # println(g_exp_ms_h)
    I2a1_vals[i] = chebyshev_shifted(exp_ms[i],k)/expm1(-2bigpi*im*g_exp_ms_h[i])
  # println(chebyshev_shifted.(exp_ms,k))
  # println("Uncorrected I2a1_vals=",I2a1_vals)
end

  debug && println("Uncorrected I2a1=",dot(fejerweights(N1,T),-exp_ms.*I2a1_vals)*im*θ/2)

  mn =  mminus+1:mplus;
  correction_coefs = Array{T}(undef,mplus-mminus)
  for (i,mni) in enumerate(mn)
    exp_msni = ap.γ(mni+hx)
    correction_coefs[i] = chebyshev_shifted(exp_msni,k)/ap.dg(exp_msni)
  end

  I2a1corr = zeros(Complex{T},N1)
  for (i,mni) in enumerate(mn)
    for j in eachindex(I2a1_vals)
      I2a1corr[j] += correction_coefs[i] / (-2bigpi*im*(g_exp_ms_h[j]-mni))
    end
  end

  fw1 = fejerweights(N1,T)
  I2a1 = zero(Complex{T})
  for i in eachindex(I2a1_vals)
    I2a1_vals[i] -= I2a1corr[i] * ap.dg(exp_ms[i])
    I2a1 += fw1[i]*(-exp_ms[i])*I2a1_vals[i]
  end
  I2a1 = I2a1*im*θ/2 + θ/2*erroradjustment(Complex{T},fejererror(N1,M1,ρ1))
  debug && println("Corrected I2a1=",I2a1)

  g_exp_mszero_h = ap.g(exp(-szero))-hx; g_exp_mszeroθ_h = ap.g(exp(-szero-im*θ))-hx
  for i in eachindex(mn)
    I2a1 += correction_coefs[i] * log((g_exp_mszeroθ_h-mn[i])/(g_exp_mszero_h-mn[i]))/(-2bigpi*im)
  end
  debug && println("Recorrected I2a1=",I2a1)


  # C2 INTEGRAL
  if nDx > 0
    I2a2 = zero(Complex{T})
    fw2 = fejerweights(N2,T)
    exp_ms = exp(-szero-Dx*(1-chebyshevpoints(N2,T)))*exp(-im*θ)
    exp_ms_sqrt = exp(-szero/2-(Dx/2)*(1-chebyshevpoints(N2,T)))*exp(-im*θ/2)
    exp_m2Dx = exp(-2Dx); exp_mDx = exp(-Dx)
    for nd = 1:nDx
      dt = zero(Complex{T})
      for i = 1:N2
        expm1_gh_invi = expm1((-2bigpi*im)*(ap.g(exp_ms[i]) - (hx + mzero))) # BIG BOPPER (~30%)
        isinf(expm1_gh_invi) ||
          (dt += fw2[i]*(-exp_ms[i])*chebyshev_shifted_sq(exp_ms_sqrt[i],k)/expm1_gh_invi)
                # BASICALLY EVERYTHING HAPPENS IN THIS CHEBYSHEV CALL (~60%)
      end
      I2a2 += Dx*(dt + erroradjustment(Complex{T},fejererror(N2,M2,ρ2)))
      exp_ms *= exp_m2Dx; exp_ms_sqrt *= exp_mDx
    end
  else
    I2a2 = zero(Complex{T})
  end
  debug && println((I2a1,I2a2))

  # TODO: BETTER XSTAR ERROR
  I2a2 += erroradjustment(Complex{T},I2a2err)

  I2 = 2real(I2a1)+2real(I2a2)

  # I1
  exp_mszero = exp(-szero)
  if k >= 2
    I1 = ((-1)^k/convert(T,k^2-1) - chebyshev_shifted(exp_mszero,k+1)/2(k+1)+chebyshev_shifted(exp_mszero,k-1)/2(k-1))/2
  elseif k == 1
    I1 = exp_mszero - exp_mszero^2
  else
    I1 = -exp_mszero
  end


  Igen = zero(T) #ap.dh(x)*(I1+I2)

  for n = ap.nmin:mzero
    γn = ap.γ(n+hx)
    Igen += chebyshev_shifted(γn,k)/ap.dg(γn)
  end

  debug &&  println(ap.dh(x)*[I1,2real(I2a1),2real(I2a2),Igen])
  # println(2real(I2a2)/(-2zeta(3)+zeta(2)- (I1+Igen+2real(I2a1))))
  -abs(ap.dh(x))*(I1+I2+Igen)
end

function transfer_basisfun_multiel{T}(ap::AbelPlanaPlan{T},x::T,K::Integer,tol=200prec(T);debug::Bool=false)
  @assert x>=0 K>=0

  α,α1,G0,G1,hx,bigpi,θ,mminus,mzero,mplus,szero,
    ρ1,M1,N1,ρ2,M2,N2,Dx,nDx,I2a2err = transfer_basisfun_parameters(ap,x,K-1,tol,debug)

  # C1 INTEGRAL
  chebpts = chebyshevpoints(N1,T)
  exp_ms = Array{Complex{T}}(undef,N1)
  I2a1_vals = Array{Complex{T}}(undef,K,N1)
  g_exp_ms_h = Array{Complex{T}}(undef,N1)
  for i = eachindex(exp_ms)
    exp_ms[i] = exp(-szero - im*θ*(1+chebpts[i])/2)
    g_exp_ms_h[i] = ap.g(exp_ms[i]) - hx
    I2a1_vals[:,i] = chebyshev_shifted_multiel(exp_ms[i],K)/expm1(-2bigpi*im*g_exp_ms_h[i])
end

  debug && println("Uncorrected I2a1=",fejerweights(N1,T)'*(-exp_ms.*transpose(I2a1_vals))*im*θ/2)

  mn =  mminus+1:mplus;
  correction_coefs = Array{T}(undef,K,mplus-mminus)
  for (i,mni) in enumerate(mn)
    exp_msni = ap.γ(mni+hx)
    correction_coefs[:,i] = chebyshev_shifted_multiel(exp_msni,K)/ap.dg(exp_msni)
  end

  I2a1corr = zeros(Complex{T},K,N1)
  for (i,mni) in enumerate(mn)
    for j = 1:N1
      I2a1corr[:,j] += correction_coefs[:,i] / (-2bigpi*im*(g_exp_ms_h[j]-mni))
    end
  end

  fw1 = fejerweights(N1,T)
  I2a1 = zeros(Complex{T},K)
  for i = 1:N1
    I2a1_vals[:,i] -= I2a1corr[:,i] * ap.dg(exp_ms[i])
    I2a1 += (fw1[i]*(-exp_ms[i]))*I2a1_vals[:,i]
  end
  I2a1 *= im*θ/2
  I2a1 += θ/2*erroradjustment(Complex{T},fejererror(N1,M1,ρ1))
  debug && println("Corrected I2a1=",I2a1)

  g_exp_mszero_h = ap.g(exp(-szero))-hx; g_exp_mszeroθ_h = ap.g(exp(-szero-im*θ))-hx
  for i in eachindex(mn)
    I2a1 += correction_coefs[:,i] * log((g_exp_mszeroθ_h-mn[i])/(g_exp_mszero_h-mn[i]))/(-2bigpi*im)
  end
  debug && println("Recorrected I2a1=",I2a1)


  # C2 INTEGRAL
  if nDx > 0
    I2a2 = zeros(Complex{T},K)
    fw2 = fejerweights(N2,T)
    exp_ms = exp(-szero-Dx*(1-chebyshevpoints(N2,T)))*exp(-im*θ)
    exp_m2Dx = exp(-2Dx)
    for nd = 1:nDx
      dt = zeros(Complex{T},K)
      for i = 1:N2
        expm1_gh_invi = expm1((-2bigpi*im)*(ap.g(exp_ms[i]) - (hx + mzero))) # BIG BOPPER (~30%)
        isinf(expm1_gh_invi) ||
          (dt += (fw2[i]*(-exp_ms[i])/expm1_gh_invi)*chebyshev_shifted_multiel(exp_ms[i],K))
                # BASICALLY EVERYTHING HAPPENS IN THIS CHEBYSHEV CALL (~60%)
      end
      I2a2 += Dx*(dt + erroradjustment(Complex{T},fejererror(N2,M2,ρ2)))
      exp_ms *= exp_m2Dx
    end
  else
    I2a2 = zeros(Complex{T},K)
  end
  debug && println((I2a1,I2a2))

  # TODO: BETTER XSTAR ERROR
  I2a2 += erroradjustment(Complex{T},I2a2err)

  I2 = 2real(I2a1)+2real(I2a2)

  # I1
  exp_mszero = exp(-szero)
  CS = chebyshev_shifted_multiel(exp_mszero,K+1)./(2(0:K))
  I1 = Array{T}(undef,K)
  K > 0 && (I1[1] = -exp_mszero)
  K > 1 && (I1[2] = exp_mszero - exp_mszero^2)
  for k = 2:K-1
    I1[k+1] = ((-1)^k/convert(T,k^2-1) - CS[k+2] + CS[k])/2
  end

  Igen = zeros(T,K) #ap.dh(x)*(I1+I2)

  for n = ap.nmin:mzero
    γn = ap.γ(n+hx)
    Igen += chebyshev_shifted_multiel(γn,K)/ap.dg(γn)
  end

  debug &&  println(ap.dh(x)*[I1 2real(I2a1) 2real(I2a2) Igen])
  # println(2real(I2a2)/(-2zeta(3)+zeta(2)- (I1+Igen+2real(I2a1))))
  -abs(ap.dh(x))*(I1+I2+Igen)
end
