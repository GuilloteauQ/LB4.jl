include("../src/LoadBalancing.jl")

function hash(s::String)::Int64
  sum(collect(s) .|> x -> Int(x))
end

function karp_rabin(scheme::Symbol, text::String, patterns::Array{String})
  n = length(patterns[1])
  m = length(text)
  hash_set = Set()
  for pattern in patterns
    push!(hash_set, hash(pattern))
  end
  count = Threads.Atomic{Int}(0)
  @LoadBalancing.lbthreads scheme for i in 1:(m-n)
    hash_sub_text = hash(text[i:(i+n-1)])
    if hash_sub_text in hash_set
      if text[i:(i+n-1)] in patterns
        Threads.atomic_add!(count, 1)
      end
    end
  end
  println(count[])
end

println("KarpRabin --------------------")
# schemes = [ :static, :gss, :tss, :dynamic, :fac2 ]
schemes = [ :static ]
for scheme in schemes
  println(scheme)
  # text = read("shakespeare.txt", String)
  # @time karp_rabin(scheme, text, ["camillo", "paulina", "leontes", "perdita", "CAMILLO", "PAULINA", "LEONTES", "PERDITA", "Camillo", "Paulina", "Leontes", "Perdita",])
  text = read("dna.200MB", String)
  @time karp_rabin(scheme, text, ["GGAT"])
end
