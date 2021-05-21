# ROW REDUCTION
function rowreduce_M{T}(M::Array{T},C,inds,fetched,qready)
  N = size(M,1)
  for i = 1:N
    i_ind = findfirst(inds,i)
    if i_ind !=0 # M ones
      v = M[:,i_ind]
      put!(C,v)
      M[i,i_ind] = 1
      for j = i+1:N
        M[j,i_ind] = 0
      end
      i_ind_p = i_ind + 1
    else
      v = fetch(C)::Vector{T}
      i_ind_p = fld(i-first(inds),step(inds))+2
    end
    put!(fetched,i)

    for k_ind = i_ind_p:length(inds) # M loop
      scale = M[i,k_ind]/v[i]
      M[i,k_ind] = scale
      for j = i+1:N
        M[j,k_ind] -= scale*v[j]
      end
    end

    take!(qready)
  end
  M
end

function rowreduce_S{T}(S::Array{T},C,inds,fetched,qready)
  N = size(S,1)
  for i = 1:N
    v = fetch(C)::Vector{T}
    put!(fetched,i)
    for k_ind in 1:length(inds) # S loop
      scale = S[i,k_ind]/v[i]
      S[i,k_ind] = scale
      for j = i+1:N
        S[j,k_ind] -= scale*v[j]
      end
    end
    take!(qready)
  end
  S
end

# BACK SUBSTITUTION
function backsub_S{T}(S::Array{T},C,inds,fetched,qready)
  # S = fetch(S_r)
  N = size(S,1)

  for i = N:-1:1
    v = fetch(C)::Vector{T}
    put!(fetched,i)
    for k_ind in 1:length(inds) # S loop
      scale = S[i,k_ind]
      for j = 1:i-1
        S[j,k_ind] -= scale*v[j]
      end
    end
    take!(qready)
  end
  S
end

function generateremoterefs(workerprocs)
  ar = Array(typeof(RemoteRef()),size(workerprocs))
  for (i,w) in enumerate(workerprocs)
    ar[i] = RemoteRef(w)
  end
  ar
end

function parallel_invert(M,workerprocs=workers();verbose=true)
  @assert size(M,1)==size(M,2)
  N = size(M,1)
  S = eye(eltype(M),N)
  nws = length(workerprocs)
  enum_workers = enumerate(workerprocs)


  verbose && println("Channels...")

  @compat C = RemoteRef(myid())
  @compat fetched_M = generateremoterefs(workerprocs)
  @compat fetched_S = generateremoterefs(workerprocs)
  @compat qready_M = generateremoterefs(workerprocs)
  @compat qready_S = generateremoterefs(workerprocs)

  # ROW REDUCTION
  verbose && println("Reducing rows...")
  @sync begin
    @async begin
      for i = 1:N
        @sync begin
          for pn = 1:nws
            @async begin
              take!(fetched_M[pn])
            end
            @async begin
              take!(fetched_S[pn])
            end
          end
        end
        take!(C)
        for pn = 1:nws
          put!(qready_M[pn],i)
          put!(qready_S[pn],i)
        end
      end
    end
    for (pn,pid) = enum_workers
      @async M[:,pn:nws:N] = remotecall_fetch(rowreduce_M,pid,M[:,pn:nws:N],C,pn:nws:N,fetched_M[pn],qready_M[pn])
      @async S[:,pn:nws:N] = remotecall_fetch(rowreduce_S,pid,S[:,pn:nws:N],C,pn:nws:N,fetched_S[pn],qready_S[pn])
    end
  end

  for pn = 1:nws
    close(fetched_M[pn])
    close(qready_M[pn])
  end

  verbose && println("Onto back substitution...")
  # BACK SUBSTITUTION
  @sync begin
    @async begin
      for i = N:-1:1
        put!(C,M[:,i])
        @sync begin
          for pn = 1:nws
            @async begin
              take!(fetched_S[pn])
            end
          end
        end
        take!(C)
        for pn = 1:nws
          put!(qready_S[pn],i)
        end
      end
    end
    for (pn,pid) = enum_workers
      @async S[:,pn:nws:N] = remotecall_fetch(backsub_S,pid,S[:,pn:nws:N],C,pn:nws:N,fetched_S[pn],qready_S[pn])
    end

  end
  close(C)
  for pn = 1:nws
    close(fetched_S[pn])
    close(qready_S[pn])
  end

  S
end
