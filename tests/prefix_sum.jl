include("../src/LB4.jl")


function prefix_sum(scheme::Symbol)
  n = 10000000
  elems = collect(1:n)

  for i in 0:(floor(Int, log2(n)) - 1)
    tmp_elems = Array{Int}(undef, n)
    @LoadBalancing.lbthreads scheme for j in 0:(n - 1)
      tmp_elems[j + 1] = elems[j + 1]
      if j >= 2^i
        tmp_elems[j + 1] += elems[j + 1 - 2^i]
      end
    end
    elems = tmp_elems
    #println(elems)
  end
end


# schemes = [ :static, :gss, :tss, :dynamic, :fac2 ]
schemes = [ :fac2 ]
for scheme in schemes
  println(scheme)
  @time prefix_sum(scheme)
end
