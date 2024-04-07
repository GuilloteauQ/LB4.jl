include("../src/LB4.jl")

function par_compute_pi(nb_iteration)::Float64
  counts = Array{Bool}(undef, nb_iteration)
  @mythreads :fac2 for i in 1:nb_iteration
    (x, y) = random_point()
    counts[i] = (x * x + y * y < 1)
  end
  4 * sum(counts) / nb_iteration
end

n = 10000000
println("PI --------------------")
print("Static ")
@time x = par_compute_pi(n)
println(x)
