# transferfun as per usual BUT also with the weight (n-m0)
function transfer_basisfun_return_parameters(em::EulerMaclaurinPlan{T},x,k::Integer,tol=200prec(T);debug::Bool=false,circratio=2) where T
  # TODO: rewrite to keep Θ biggish at cost of large m̂

  # ITERATION UP TO EM REGION CONSTANTS
  L = one(T)#max(one(T),T(k)/40) #one(T) #whatever
  srd_sqrt = sqrt(em.saferad)
  Θ = min(asin(min(L/max(k,1)*(cos(em.Θmax)/srd_sqrt-srd_sqrt),1)),
        em.Θmax)

  # Θ = em.Θmax/max(k,1) ## TODO: better?
  debug && println("L=",L," Θ=",Θ)

  rmax = min(em.saferad,
    em.α1 == 0 ? (em.G0/em.G1 * sin(Θ/2))^(1/em.α) : ((abs(em.α1)/em.α)*em.G0/em.G1 * sin(Θ/2))^(1/(em.α-em.α1)) )
  rdmin = g0int(em,rmax) + g1int(em,rmax)
  Θd = min(em.α*Θ - asin(g1int(em,rmax)/g0int(em,rmax)), T(3)/2) #slightly smaller than pi/2
  debug && println("rmax=",rmax," rdmin=",rdmin," Θd=",Θd)

  # CIRCLE CONSTANTS
  p = roundupbound(-log(tol)/2) #TODO: is this /2 thing part of circratio?
  rd_m0 = 2p/T(pi)
  debug && println("p=",p," rd_m0=",rd_m0," circratio=",circratio)

  # ITERATION UP TO EM REGION CONSTANTS
  m_orig = em.h(x)+T(em.nmin)
  m̂_discrim = (2T(π) + rdmin)^2-imag(m_orig)^2
  m̂ = m̂_discrim ≤ 0 ? 0 : roundupbound(max(0,sqrt(m̂_discrim)-real(m_orig)))
  m̂ = max(m̂,roundupbound(max(0,cot(Θd)*abs(imag(m_orig))-real(m_orig)+rd_m0*csc(Θd))))
  m0 = m_orig + m̂
  debug && println(" m_orig=",m_orig," m̂=",m̂," m0=",m0)

  # ERROR BOUND CONSTANTS
  ϕmax = exp(L)*rd_m0
  ϕbd = exp(L)/(em.G0 * rmax^(-em.α-1) - em.G1 * rmax^(-em.α1-1))
  κ = 1
  debug && println("ϕmax=",ϕmax," ϕbd=",ϕbd," κ=",κ)

  cc = ((abs(m0)-rdmin)^(1-2p+κ)/(real(m0)/abs(m0))+(sin(Θd)*real(m0)-cos(Θd)*imag(m0))^(1-2p+κ)/sin(Θd))/(2p-1-κ)
  debug && println("cc=",cc)

  # CIRCLE TAYLOR TRANSFORMATION CONSTANTS
  Ncirc = nextpow(2,max(2p,roundupbound(log(max(1,ϕmax)/tol + 1)/log(circratio))))
  Mcirc = ϕmax #TODO
  debug && println("Ncirc=",Ncirc," Mcirc=",Mcirc)

  return p, Ncirc, Mcirc, rd_m0,circratio,ϕmax,ϕbd,m0,cc,m_orig,m̂
end

# does transferfun as per usual BUT also with the weight (n-m0)
function transfer_basisfun_withreturn_multiel_em(ir::IntermittentReturnMap{T},x::Complex{T},K::Integer,
      tol,params_em,params_emr;debug::Bool=false) where T
  em = ir.em
  (p, Ncirc, Mcirc, rd_m0,circratio,ϕmax,ϕbd,m0,cc) = params_em
  (_p, Ncircr, Mcircr, _rd_m0,_circratio,ϕmaxr,ϕbdr,_m0,ccr) = params_emr
  Ncirc = max(Ncirc,Ncircr)

  ϕoncircle = Array{Complex{T}}(undef,K,Ncirc)
  # ϕroncircle = Array{Complex{T}}(undef,K,Ncirc)
  z = complex(rd_m0/circratio); twid = exp(2im*T(pi)/Ncirc)
  for i = 1:Ncirc
    γv = em.γ(m0+z)
    ϕoncircle[:,i] = chebyshev_shifted_multiel(γv,K)/em.dg(γv)
    # ϕroncircle[:,i] = z*chebyshev_shifted_multiel(γv,K)/em.dg(γv)
    z *= twid
  end
  Taylorϕ = Array{Complex{T}}(undef,p,K)
  Taylorϕr = Array{Complex{T}}(undef,p,K)
  for k = 1:K
    tt= taylortransform(vec(ϕoncircle[k,:]),rd_m0/circratio,Mcirc,rd_m0)
    # ttr= taylortransform(vec(ϕroncircle[k,:]),rd_m0/circratio,Mcirc,rd_m0)
    Taylorϕ[:,k] = tt[2:2:2p]
    Taylorϕr[:,k] = tt[1:2:2p-1]#ttr[2:2:2p]
  end
  B2k = bernoulli2k(p,T)
  Icirc = zeros(Complex{T},K)
  Icircr = zeros(Complex{T},K)
  for k = 1:K
    for i = 1:p
      Icirc[k] -= B2k[i+1]*Taylorϕ[i,k]/(2i)
      Icircr[k] -= B2k[i+1]*Taylorϕr[i,k]/(2i)
    end
  end
  debug && println("Icirc=",Icirc," Icircr=",Icircr)

  # ERROR
  Ierr= erroradjustment(Complex{T},abs(B2k[end]) * ϕbd * cc)
  Ierr_r= erroradjustment(Complex{T},abs(B2k[end]) * ϕbdr * ccr)
  debug && println("Ierr=",Ierr," Ierr_r=",Ierr_r)

  # HALF-COEFFICIENT SUM
  γm0 = em.γ(m0)
  Ihconstantcoeff = chebyshev_shifted_multiel(γm0,K)#
  for i = 1:K # weird bug
    Ihconstantcoeff[i] /= 2em.dg(γm0)
  end
  Ihconstantcoeffr = zero(Complex{T})
  debug && println("Half-coefficient=",Ihconstantcoeff,"Half-coefficient for return=",Ihconstantcoeffr)

  # NORMAL INTEGRAL
  CS = chebyshev_shifted_multiel(γm0,K+1)./(2(0:K))
  Iintegral = Array{Complex{T}}(undef,K)
  K > 0 && (Iintegral[1] = -γm0)
  K > 1 && (Iintegral[2] = γm0 - γm0^2)
  for k = 2:K-1
    Iintegral[k+1] = ((-1)^k/convert(T,k^2-1) - CS[k+2] + CS[k])/2
  end
  debug && println("Iintegral=",Iintegral)

  # THE OTHER ONE
  Iintegralr = returnintegral_multiel(ir.a,complex(γm0),K)#+em.h(x)*Iintegral
  debug && println("Iintegralr=",Iintegralr)

  dhx= em.dh(x) * em.σh
  debug && println("σh*dh(x)=",dhx)
  return -(Icirc + Ihconstantcoeff + Iintegral + Ierr) * dhx,
    -(Icircr + Ihconstantcoeffr + Iintegralr + Ierr_r) * dhx
end
function transfer_basisfun_withreturn_multiel_em(ir::IntermittentReturnMap{T},x::T,args...;kwargs...) where T<:Real
  a,b=transfer_basisfun_withreturn_multiel_em(ir,complex(x),args...;kwargs...)
  real(a), real(b)
end

function returnintegral_multiel(a::AbelFunction{T},x::Complex{T},K) where T
  NT = zeros(Complex{T},K)
  x̂ = x^(a.r.α)
  logx = log(x)
  gcoeffs0 = a.coeffs[1]
  gcoeffsl = a.coeffs[2]*a.r.α
  gcoeffs = [a.offset-a(x);a.coeffs[3:end]]
  for k = 0:K-1
    # chebcoefs = Array{T}(undef,k+1)
    chebcoef_lx = (-1)^k * x
    for l = 0:k
      gsum =  gcoeffs0 / x̂ / (l+1-a.r.α)
      gsum+=  gcoeffsl / (l+1)^2 * ((l+1)*logx-1)
      x̂pow = one(Complex{T})
      for n = 0:length(gcoeffs)-1
        gsum += gcoeffs[n+1] * x̂pow / (l+1+a.r.α*n)
        x̂pow *= x̂
      end
      NT[k+1] += gsum*chebcoef_lx
      chebcoef_lx *= -(x*(4(k+l)*(k-l)))/((2l+1)*(2l+2))
    end
  end
  -NT
end

function _transferloop_withreturn_multiel!(r,K,x̂n,vn,mord,transfermass,transfermassr)
  xnm = x̂n^(1/r.α)#AbelFunctions.unhat(x̂n,r)
  tm = chebyshev_shifted_multiel(xnm*r.sgn+r.p,K) *(xnm/x̂n * vn)
  @inbounds transfermass[:] += tm
  @inbounds transfermassr[:] += mord * tm
  x̂n = AbelFunctions.mapinv_trans(r,x̂n)
  vn /= AbelFunctions.mapD_trans(r,x̂n)
  # println(x̂n)
  x̂n, vn, mord + 1
end


function transfer_basisfun_withreturn_multiel(irm::IntermittentReturnMap{T},x,K::Integer,tol=200prec(T);debug::Bool=false,circratio=2) where T
  transfermass = zeros(T,K)
  transfermassr = zeros(T,K)

  mord = 0
  x̂n = AbelFunctions.hat(x,irm.r)
  vn = x̂n/(irm.r.sgn*(x-irm.r.p))
  ŝr̂ = irm.em.saferad^(irm.r.α)
  while ~(abs(x̂n) < ŝr̂)
    x̂n,vn,mord = _transferloop_withreturn_multiel!(irm.r,K,x̂n,vn,mord,transfermass,transfermassr)
  end
  debug && println("number of iterations to AbelFunction noptrad = $mord, transfermass = $transfermass, transfermassr = $transfermassr")
  # with diameter $(broadcast(diam,transfermass))")
  xnm_em = x̂n^(1/irm.r.α)
  xn_em = irm.r.p + irm.r.sgn*xnm_em#AbelFunctions.unhat(x̂n,irm.r)
  params = transfer_basisfun_parameters(irm.em,xn_em,K-1,tol;debug=debug,circratio=circratio)
  paramsr = transfer_basisfun_return_parameters(irm.em,xn_em,K-1,tol;debug=debug,circratio=circratio)
  m̂ = params[end]

  debug && println("E-M m̂ = $(params[end])")
  tm_em,tm_emr = transfer_basisfun_withreturn_multiel_em(irm,xn_em,K,tol,params[1:end-2],paramsr[1:end-2],debug=debug)
  tm_em *= (vn * xnm_em/x̂n);  tm_emr *= (vn * xnm_em/x̂n)
  tm_emr += (m̂+mord)*tm_em
  transfermass += tm_em; transfermassr += tm_emr;
  debug && println("E-M transfer mass = $tm_em, transfer mass with return weight = $tm_emr")#, with diameter $(broadcast(diam,tm_em))")
  for i = 1:m̂
    x̂n,vn,mord = _transferloop_withreturn_multiel!(irm.r,K,x̂n,vn,mord,transfermass,transfermassr)
    # debug && mord += 1
  end

  debug && println("number of iterations to E-M point = $(mord), ",
    "total transfermass = $transfermass, total transfermassr = $transfermassr")
  debug && println("E-M point at $(AbelFunctions.unhat(x̂n,irm.r))")
  # debug && println("total transfermass = $transfermass")

  transfermassr += transfermass # since return time starts at 0 rather than 1
  transfermass, transfermassr
end
