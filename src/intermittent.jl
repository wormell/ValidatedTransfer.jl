@compat struct AbelDerivative{T,ffa,dffa}; a::AbelFunction{T,ffa,dffa}; end
@compat (dg::AbelDerivative)(x) = mapD(dg.a,x)
@compat struct AbelInverse{T,ffa,dffa}; a::AbelFunction{T,ffa,dffa}; end
@compat (γ::AbelInverse)(x) = mapinv(γ.a,x)

@compat struct IntermittentReturnMap{T<:Real,ffa,dffa} # full domain = [0,1] as per usual
  r::NeutralRecurrence{T,ffa,dffa}
  a::AbelFunction{T,ffa,dffa}
  em::EulerMaclaurinPlan{T,AbelFunction{T,ffa,dffa},AbelDerivative{T,ffa,dffa},AbelInverse{T,ffa,dffa},
    AbelFunction{T,ffa,dffa},AbelDerivative{T,ffa,dffa}} # eep
end
function IntermittentReturnMap(r::NeutralRecurrence{T,ffa,dffa},#returnset::Interval{T}=Interval(T(1/2),T(1)),
    a=AbelFunction(r,basepoint=T(Int(r.sgn>0)))) where {T,ffa,dffa}
  @assert (r.p == 0 && r.sgn == 1) || (r.p == 1 && r.sgn == -1) # so it's Markov on [0,1]
  em=generateemplan(a)
  IntermittentReturnMap{T,ffa,dffa}(r,a,em)
end
IntermittentReturnMap(a::AbelFunction{T,ffa,dffa}) where {T,ffa,dffa} = IntermittentReturnMap(a.r,a)

function generateemplan(a::AbelFunction{T}) where T
  # determining G1 vals: TODO check this
  G1 = T(0); G1d = abs(a.coeffs[2])
  crt = one(T)
  for i = 1:a.ncall
    crt *= a.noptrad
    crtp = crt * abs(a.coeffs[i])
    G1 += crtp
    G1d += i*crtp
  end
  saferad = (3a.noptrad/4)^(1/a.r.α) #3/4 just some number
  Θmax = T(π)/4a.r.α
  badcircpt = (saferad*exp(im*Θmax))^a.r.α
  G1 = max(abs(a.coeffs[2])+(G1+AbelFunctions.abelerror(a,badcircpt))/log(saferad),
    a.r.α*G1d+ AbelFunctions.abelerrorD(a,badcircpt)*abs(badcircpt)*a.r.α)
  G0 = a.coeffs[1]*a.r.α
  EulerMaclaurinPlan(a,AbelDerivative(a),AbelInverse(a),a,AbelDerivative(a),a.r.α,zero(T),G0,G1;
        nmin=0,saferad=saferad,Θmax=Θmax,σh=-1)
end


# TRANSFERRING

function transferloop(r,k,x̂n,vn,transfermass)
  xn = AbelFunctions.unhat(x̂n,r)
  dn = x̂n^(1/r.α-1)
  transfermass += chebyshev_shifted(xn,k) *(dn * vn)
  x̂n = AbelFunctions.mapinv_trans(r,x̂n)
  vn /= AbelFunctions.mapD_trans(r,x̂n)
  # println(x̂n)
  x̂n, vn, transfermass
end

function transfer_basisfun(irm::IntermittentReturnMap{T},x,k::Integer,tol=200prec(T);debug::Bool=false,circratio=2) where T
  transfermass = zero(typeof(x))

  mord = 0
  x̂n = AbelFunctions.hat(x,irm.r)
  vn = x̂n/(irm.r.sgn*(x-irm.r.p))
  ŝr̂ = irm.em.saferad^(irm.r.α)
  while ~(abs(x̂n) < ŝr̂)
    x̂n,vn,transfermass = transferloop(irm.r,k,x̂n,vn,transfermass)
    mord += 1
  end
  debug && println("number of iterations to AbelFunction noptrad = $mord, transfermass = $transfermass")
  # with diameter $(broadcast(diam,transfermass))")
  xnm_em = x̂n^(1/irm.r.α)
  xn_em = irm.r.p + irm.r.sgn*xnm_em#AbelFunctions.unhat(x̂n,irm.r)
  params = transfer_basisfun_parameters(irm.em,xn_em,k,tol;debug=debug,circratio=circratio)
  debug && println("E-M m̂ = $(params[end])")
  tm_em = transfer_basisfun_em(irm.em,xn_em,k,tol,params[1:end-2],debug=debug) * (vn * xnm_em/x̂n)
  transfermass += tm_em
  debug && println("E-M transfer mass = $tm_em")#, with diameter $(broadcast(diam,tm_em))")
  for i = 1:params[end] # = m̂
    x̂n,vn,transfermass = transferloop(irm.r,k,x̂n,vn,transfermass)
    # debug && mord += 1
  end

  debug && println("number of iterations to E-M point = $(mord+params[end]), total transfermass = $transfermass")
  debug && println("E-M point at $(AbelFunctions.unhat(x̂n,irm.r))")
  # debug && println("total transfermass = $transfermass")

  transfermass
end

function _transferloop_multiel!(r,K,x̂n,vn,transfermass)
  xnm = x̂n^(1/r.α)#AbelFunctions.unhat(x̂n,r)
  @inbounds transfermass[:] += chebyshev_shifted_multiel(xnm*r.sgn+r.p,K) *(xnm/x̂n * vn)
  x̂n = AbelFunctions.mapinv_trans(r,x̂n)
  vn /= AbelFunctions.mapD_trans(r,x̂n)
  # println(x̂n)
  x̂n, vn
end


function transfer_basisfun_multiel(irm::IntermittentReturnMap{T},x,K::Integer,tol=200prec(T);debug::Bool=false,circratio=2) where T
  transfermass = zeros(T,K)

  mord = 0
  x̂n = AbelFunctions.hat(x,irm.r)
  vn = x̂n/(irm.r.sgn*(x-irm.r.p))
  ŝr̂ = irm.em.saferad^(irm.r.α)
  while ~(abs(x̂n) < ŝr̂)
    x̂n,vn = _transferloop_multiel!(irm.r,K,x̂n,vn,transfermass)
    debug && (mord += 1)
  end
  debug && println("number of iterations to AbelFunction noptrad = $mord, transfermass = $transfermass")
  # with diameter $(broadcast(diam,transfermass))")
  xnm_em = x̂n^(1/irm.r.α)
  xn_em = irm.r.p + irm.r.sgn*xnm_em#AbelFunctions.unhat(x̂n,irm.r)
  params = transfer_basisfun_parameters(irm.em,xn_em,K-1,tol;debug=debug,circratio=circratio)
  debug && println("E-M m̂ = $(params[end])")
  tm_em = transfer_basisfun_multiel_em(irm.em,xn_em,K,tol,params[1:end-2],debug=debug) * (vn * xnm_em/x̂n)
  transfermass += tm_em
  debug && println("E-M transfer mass = $tm_em")#, with diameter $(broadcast(diam,tm_em))")
  for i = 1:params[end] # = m̂
    x̂n,vn = _transferloop_multiel!(irm.r,K,x̂n,vn,transfermass)
    # debug && mord += 1
  end

  debug && println("number of iterations to E-M point = $(mord+params[end]), total transfermass = $transfermass")
  debug && println("E-M point at $(AbelFunctions.unhat(x̂n,irm.r))")
  # debug && println("total transfermass = $transfermass")

  transfermass
end
