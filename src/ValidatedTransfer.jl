module ValidatedTransfer
using IntervalArithmetic, Compat, FastTransforms, JLD, FileIO

include("general.jl")
include("transforms.jl")
include("transferbounds.jl")

function filltransfer(J,K,transferfn_N,helperdata,aliaserror,entrybound,N_interp=maximum(J),T=BigInterval)
  S = Array{T}(length(J),length(K))
  for (k_ind,k) in enumerate(K)
    B = Array{T}(N_interp)
    for i in eachindex(B)
        B[i] = transferfn_N(i,k,helperdata)
    end
    S[:,k_ind] = chebyshevtransform(B)[J]
    for j in J
      ae = aliaserror(j,k,N_interp)
      S[j,k_ind] += erroradjustment(Interval(ae))
      S[j,k_ind] = S[j,k_ind] ∩ erroradjustment(entrybound(j,k))
    end
  end

  S
end

function maketransfer(N,transferfn,helperdata,aliaserror,entrybound,workerprocs=workers(),N_interp=N,T=BigInterval)
  nws = length(workerprocs)
  enum_workers = enumerate(workerprocs)

  M = Array(T,N,N)
  @sync begin
    for (pn,pid) in enum_workers
      @async M[:,pn:nws:N] = remotecall_fetch(filltransfer,pid,1:N,pn:nws:N,transferfn,helperdata,aliaserror,entrybound,N_interp,T)
    end
  end
  M
end

function makesolutioninv!(M::Array)
  # M = fetch(M_r)
  for k in 1:size(M,2)
    for j = 1:size(M,1)
      M[j,k] *= -1
    end

    mod(k,2) == 1 && (k== 1 ? (M[1,k] += 1) :
        (M[1,k] -= 1/convert(eltype(M),(k-1)^2-1)) )
    M[k,k] += 1
  end
  M
end

include("inversion.jl")

function getsolution!(M,workerprocs=workers();verbose=true)
  makesolutioninv!(M)
  parallel_invert(M,workerprocs;verbose=verbose)
end

function getsolution(M,workerprocs=workers();verbose=true)
  K=makesolutioninv!(deepcopy(M))
  parallel_invert(K,workerprocs;verbose=verbose)
end

include("files.jl")

export BigInterval, erroradjustment
export ChebyshevZeroOneBV, AnalyticEntryBounds, FirstOrderLYBounds
export logSEnorm, logEnorm, Snorm, aliaserror, entrybound
export transfer_basisfun, transfer_basisfun_multiel
export maketransfer, getsolution
export savezip, loadzip

end # module
