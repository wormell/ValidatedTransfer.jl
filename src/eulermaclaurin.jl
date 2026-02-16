struct EulerMaclaurinPlan{T,Tg,Tdg,Tγ,Th,Tdh}
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
  σh::Int
end
EulerMaclaurinPlan(args...;kwargs...) = _EulerMaclaurinPlan(args...;kwargs...)

function _EulerMaclaurinPlan(g,dg,γ,h,dh,α,α1,G0,G1;nmin::Int=0,saferad=1,Θmax=2,σh=intervalsign(dh(saferad/2)))
  (α,α1,G0,G1,saferad,Θmax) = promote(α,α1,G0,G1,saferad,Θmax)
  T = typeof(α)
  G1 = max(G1,G0*prec(T))
   Θmax = min(Θmax,T(pi)/2α,T(pi)/3)
  EulerMaclaurinPlan{T,typeof(g),typeof(dg),typeof(γ),typeof(h),typeof(dh)}(g,dg,γ,h,dh,α,α1,G0,G1,nmin,saferad,Θmax,σh)
end


g0int(em::EulerMaclaurinPlan,rmax) = em.G0 * rmax^(-em.α) / em.α
g1int(em::EulerMaclaurinPlan,rmax) = em.G1 * ( em.α1==0 ? log(rmax) : rmax^(-em.α1)/abs(em.α1) )
function transfer_basisfun_parameters(em::EulerMaclaurinPlan{T},x,k::Integer,tol=200prec(T);debug::Bool=false,circratio=2) where T
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
  ϕmax = exp(L)
  ϕbd = ϕmax/(em.G0 * rmax^(-em.α-1) - em.G1 * rmax^(-em.α1-1))
  debug && println("ϕmax=",ϕmax," ϕbd=",ϕbd)

  cc = ((abs(m0)-rdmin)^(1-2p)/(real(m0)/abs(m0))+(sin(Θd)*real(m0)-cos(Θd)*imag(m0))^(1-2p)/sin(Θd))/(2p-1)
  debug && println("cc=",cc)

  # CIRCLE TAYLOR TRANSFORMATION CONSTANTS
  Ncirc = nextpow(2,max(2p,roundupbound(log(max(1,ϕmax)/tol + 1)/log(circratio))))
  Mcirc = ϕmax #TODO
  debug && println("Ncirc=",Ncirc," Mcirc=",Mcirc)

  return p, Ncirc, Mcirc, rd_m0,circratio,ϕmax,ϕbd,m0,cc,m_orig,m̂
end

function transfer_basisfun_em(em::EulerMaclaurinPlan{T},x::Complex{T},k::Integer,tol,params_em;debug::Bool=false) where T
  (p, Ncirc, Mcirc, rd_m0,circratio,ϕmax,ϕbd,m0,cc) = params_em
  ϕoncircle = Array{Complex{T}}(undef,Ncirc)
  z = complex(rd_m0/circratio); twid = exp(2im*T(pi)/Ncirc)
  for i = 1:Ncirc
    γv = em.γ(m0+z)
    ϕoncircle[i] = chebyshev_shifted(γv,k)/em.dg(γv)
    z *= twid
  end
  Taylorϕ = taylortransform(ϕoncircle,rd_m0/circratio,Mcirc,rd_m0)
  println(Taylorϕ[1:10])
  B2k = bernoulli2k(p,T)
  Icirc = zero(Complex{T})
  for i = 1:p
    Icirc -= B2k[i+1]*Taylorϕ[2i]/(2i)
  end
  debug && println("Icirc=",Icirc)

  # ERROR
  Ierr= erroradjustment(Complex{T},abs(B2k[end]) * ϕbd * cc)
  debug && println("Ierr=",Ierr)

  # INTEGRAL
  γm0 = em.γ(m0)
  if k >= 2
    Iintegral = ((-1)^k/convert(T,k^2-1) - chebyshev_shifted(γm0,k+1)/2(k+1)+chebyshev_shifted(γm0,k-1)/2(k-1))/2
  elseif k == 1
    Iintegral = γm0 - γm0^2
  elseif k == 0
    Iintegral = -γm0
  end
  debug && println("Iintegral=",Iintegral)

  # HALF CONSTANT
  Ihconstantcoeff = chebyshev_shifted(γm0,k)/em.dg(γm0)/2
  debug && println("Half-coefficient=",Ihconstantcoeff)

  return -(Icirc + Ierr + Ihconstantcoeff + Iintegral)*em.dh(x)*em.σh
end
transfer_basisfun_em(em::EulerMaclaurinPlan{T},x::T,args...;kwargs...) where T<:Real =
  real(transfer_basisfun_em(em,complex(x),args...;kwargs...))

function transfer_basisfun(em::EulerMaclaurinPlan{T},x::Complex{T},k::Integer,tol=200prec(T);debug::Bool=false) where T
  params = transfer_basisfun_parameters(em,x,k,tol;debug=debug)
  em_tot= transfer_basisfun_em(em,x,k,tol,params[1:end-2],debug=debug)
  m_orig,m̂ = params[end-1:end]

  # TAIL SUMS
  Iconstcoeffs = zero(Complex{T})
  for i = 0:m̂-1
    γv = em.γ(m_orig + i)
    Iconstcoeffs += chebyshev_shifted(γv,k)/em.dg(γv)
  end
  debug && println("Iconstcoeffs=",Iconstcoeffs)

  return -Iconstcoeffs*em.dh(x)*em.σh + em_tot
end


function transfer_basisfun_multiel_em(em::EulerMaclaurinPlan{T},x::Complex{T},K::Integer,tol,params_em;debug::Bool=false) where T
  p, Ncirc, Mcirc, rd_m0,circratio,ϕmax,ϕbd,m0,cc = params_em
  ϕoncircle = Array{Complex{T}}(undef,K,Ncirc)
  z = complex(rd_m0/circratio); twid = exp(2im*T(pi)/Ncirc)
  for i = 1:Ncirc
    γv = em.γ(m0+z)
    ϕoncircle[:,i] = chebyshev_shifted_multiel(γv,K)/em.dg(γv)
    z *= twid
  end
  Taylorϕ = Array{Complex{T}}(undef,p,K)
  for k = 1:K
    Taylorϕ[:,k] = taylortransform(vec(ϕoncircle[k,:]),rd_m0/circratio,Mcirc,rd_m0)[2:2:2p]
  end
  B2k = bernoulli2k(p,T)
  Icirc = zeros(Complex{T},K)
  for k = 1:K
    for i = 1:p
      Icirc[k] -= B2k[i+1]*Taylorϕ[i,k]/(2i)
    end
  end
  debug && println("Icirc=",Icirc)

  # ERROR
  Ierr= erroradjustment(Complex{T},abs(B2k[end]) * ϕbd * cc)
  debug && println("Ierr=",Ierr)

  # HALF-COEFFICIENT SUM
  γm0 = em.γ(m0)
  Ihconstantcoeff = chebyshev_shifted_multiel(γm0,K)/em.dg(γm0)/2
  debug && println("Half-coefficient=",Ihconstantcoeff)

  # INTEGRAL
  CS = chebyshev_shifted_multiel(γm0,K+1)./(2(0:K))
  Iintegral = Array{Complex{T}}(undef,K)
  K > 0 && (Iintegral[1] = -γm0)
  K > 1 && (Iintegral[2] = γm0 - γm0^2)
  for k = 2:K-1
    Iintegral[k+1] = ((-1)^k/convert(T,k^2-1) - CS[k+2] + CS[k])/2
  end
  debug && println("Iintegral=",Iintegral)
  return -(Icirc + Ihconstantcoeff + Iintegral .+ Ierr) * (em.dh(x) * em.σh)
end
transfer_basisfun_multiel_em(em::EulerMaclaurinPlan{T},x::T,args...;kwargs...) where T<:Real =
  real(transfer_basisfun_multiel_em(em,complex(x),args...;kwargs...))

function transfer_basisfun_multiel(em::EulerMaclaurinPlan{T},x,K::Integer,tol=200prec(T);debug::Bool=false) where T
  params = transfer_basisfun_parameters(em,x,K-1,tol;debug=debug)
  em_tot = transfer_basisfun_multiel_em(em,x,K,tol,params[1:end-2];debug=debug)
  m_orig,m̂ = params[end-1:end]

  # COEFFICIENT SUMS
  Iconstcoeffs = zeros(Complex{T},K)
  for i = 0:m̂-1
    γv = em.γ(m_orig + i)
    Iconstcoeffs += chebyshev_shifted_multiel(γv,K)/em.dg(γv)
  end
  debug && println("Iconstcoeffs=",Iconstcoeffs)

  return -Iconstcoeffs  * (em.dh(x) * em.σh) + em_tot
end

# GENERAL FUNS

# function transfer_fun_parameters(em::EulerMaclaurinPlan{T},x,ρ,M,tol=200prec(T);debug::Bool=false) where T
#   ϕbd =
#   srd_sqrt = sqrt(em.saferad)
#   Θ = min(asin(min(L/max(k,1)*(cos(em.Θmax)/srd_sqrt-srd_sqrt),1)),
#         em.Θmax)
#
#   # Θ = em.Θmax/max(k,1) ## TODO: better?
#   debug && println("L=",L," Θ=",Θ)
#
#   rmax = min(em.saferad,
#     em.α1 == 0 ? (em.G0/em.G1 * sin(Θ/2))^(1/em.α) : ((abs(em.α1)/em.α)*em.G0/em.G1 * sin(Θ/2))^(1/(em.α-em.α1)) )
#   rdmin = g0int(em,rmax) + g1int(em,rmax)
#   Θd = min(em.α*Θ - asin(g1int(em,rmax)/g0int(em,rmax)), T(3)/2) #slightly smaller than pi/2
#   debug && println("rmax=",rmax," rdmin=",rdmin," Θd=",Θd)
#
#   ϕmax = exp(L)
#   ϕbd = 1/(em.G0 * rmax^(-em.α-1) - em.G1 * rmax^(-em.α1-1))
#   debug && println("ϕmax=",ϕmax," ϕbd=",ϕbd)
#
#   p = roundupbound(-log(tol)/2) #TODO: is this /2 thing part of circratio?
#   rd_m0 = 2p/T(pi)
#   circratio = 2
#   Ncirc = nextpow2(max(2p,roundupbound(log(max(1,ϕmax)/tol + 1)/log(circratio))))
#   Mcirc = ϕmax #TODO
#   debug && println("p=",p," rd_m0=",rd_m0," circratio=",circratio," Ncirc=",Ncirc," Mcirc=",Mcirc)
#
#   m_orig = em.h(x)+T(em.nmin)
#   m̂_discrim = (2T(π) + rdmin)^2-imag(m_orig)^2
#   m̂ = m̂_discrim ≤ 0 ? 0 : roundupbound(max(0,sqrt(m̂_discrim)-real(m_orig)))
#   m̂ = max(m̂,roundupbound(max(0,cot(Θd)*abs(imag(m_orig))-real(m_orig)+rd_m0*csc(Θd))))
#   m0 = m_orig + m̂
#   debug && println(" m_orig=",m_orig," m̂=",m̂," m0=",m0)
#
#   cc = ((abs(m0)-rdmin)^(1-2p)/(real(m0)/abs(m0))+(sin(Θd)*real(m0)-cos(Θd)*imag(m0))^(1-2p)/sin(Θd))/(2p-1)
#   debug && println("cc=",cc)
#
#   return p, Ncirc, Mcirc, rd_m0,circratio,ϕmax,ϕbd,m0,cc,m_orig,m̂
# end
