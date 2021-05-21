function savezip{T}(fname,S::Array{IntervalArithmetic.Interval{T}},extras...)
  M_lo = Array{T}(undef,size(S))
  M_hi = Array{T}(undef,size(S))
  for ind in eachindex(M_lo)
    M_lo[ind] = S[ind].lo
    M_hi[ind] = S[ind].hi
  end
  tpath = mktempdir()*"/"*splitext(splitdir(fname)[2])[1]*".jld"
  save(File(format"JLD",tpath),"lo",M_lo,"hi",M_hi,extras...)
  run(`zip -jr9 $fname $tpath`)
end

function loadzip(fname)
  tpath = mktempdir()*"/"
  run(`unzip -j $fname -d $tpath`)

  ld = load(tpath*splitext(splitdir(fname)[2])[1]*".jld")
  run(`rm -r $tpath`)
  T = typeof(ld["lo"][1])
  J,K = size(ld["lo"])
  M = Array{IntervalArithmetic.Interval{T}}(J,K)
  for j = 1:J, k = 1:K
    M[j,k] = IntervalArithmetic.Interval(ld["lo"][j,k],ld["hi"][j,k])
  end
  M
end
