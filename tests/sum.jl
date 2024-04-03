include("../src/LoadBalancing.jl")
using BenchmarkTools

function random_point()
  (rand(), rand())
end

function par_compute_pi(nb_iteration, scheme::Symbol)::Float64
  counts = Array{Bool}(undef, nb_iteration)
  @LoadBalancing.lbthreads scheme for i in 1:nb_iteration
    (x, y) = random_point()
    counts[i] = (x * x + y * y < 1)
  end
  4 * sum(counts) / nb_iteration
end

println("PI --------------------")
n = 10000000
# schemes = [ :static, :gss, :tss, :dynamic, :fac2 ]
schemes = [:gss]
for scheme in schemes
  print(scheme)
  @time par_compute_pi(n, scheme)
end

