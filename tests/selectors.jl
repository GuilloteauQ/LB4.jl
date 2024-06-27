include("../src/LB4.jl")

function diverges(c::Complex, nb_iterations::Int64, epsilon::Float64)::Bool
  z = 0 + im * 0
  i = 0
  while abs(z) < epsilon && i < nb_iterations
    z = z * z + c
    i = i + 1
  end
  i < nb_iterations
end

function mandelbrot(scheme::Symbol)
  xmin = -4.0
  xmax =  4.0
  ymin = -1.0
  ymax =  1.0
  delta = 0.005

  nb_rows = floor(Int64, (ymax - ymin) / delta)
  nb_cols = floor(Int64, (xmax - xmin) / delta)

  data = Matrix{Bool}(undef, nb_rows, nb_cols)
  total_iters = nb_cols * nb_rows

  plop = @LB4.lbthreads scheme for k in 0:(nb_cols * nb_rows - 2)
    y = div(k, nb_cols) 
    x = (k) % nb_cols
    c = ((xmax - xmin) * x / nb_cols + xmin) + im * ((ymax - ymin) * y / nb_rows + ymin)
    data[y + 1, x + 1] = diverges(c, 10000, 2.0)
  end
  plop
end

mutable struct UCB1
  info_times::Vector{Float64}
  info_nb_choice::Vector{Int64}
  schemes::Vector{Symbol}
  last_choice::Int64
  iteration::Int64

  UCB1(schemes) = new(zeros(Float64, length(schemes)), zeros(Int64, length(schemes)), schemes, -1, 0)
end

function choice!(selector::UCB1)
  nb_times, index = findmin(selector.info_nb_choice)
  if nb_times != 0
    selector.iteration += 1
    _, index = (zip(selector.info_times, selector.info_nb_choice) .|> x -> x[1]/x[2] - sqrt(0.1 * log(selector.iteration) / x[2])) |> findmin
  end
  selector.info_nb_choice[index] += 1
  selector.last_choice = index
  return selector.schemes[index]
end

function update!(selector::UCB1, perf)
  selector.info_times[selector.last_choice] += perf
end

mutable struct UCBV
  info_times::Vector{Float64}
  info_times_sq::Vector{Float64}
  info_nb_choice::Vector{Int64}
  schemes::Vector{Symbol}
  last_choice::Int64
  iteration::Int64

  UCBV(schemes) = new(zeros(Float64, length(schemes)), zeros(Float64, length(schemes)), zeros(Int64, length(schemes)), schemes, -1, 0)
end

function choice!(selector::UCBV)
  nb_times, index = findmin(selector.info_nb_choice)
  if nb_times != 0
    selector.iteration += 1
    _, index = (zip(selector.info_times, selector.info_times_sq, selector.info_nb_choice) .|> x -> x[1]/x[3] - (sqrt(1 * log(selector.iteration) * (x[2]/x[3] - (x[1]/x[3])^2) / x[3]) + 3 * log(selector.iteration) / x[3])) |> findmin
  end
  selector.info_nb_choice[index] += 1
  selector.last_choice = index
  return selector.schemes[index]
end

function update!(selector::UCBV, perf)
  selector.info_times[selector.last_choice] += perf
  selector.info_times_sq[selector.last_choice] += perf*perf
end


function main()
  schemes = [ :static, :gss, :tss, :dynamic, :fac2, :viss, :rnd, :auto ]

  selector = UCB1(schemes)

  for t in 1:1000
    scheme = choice!(selector)
    runtime_data = mandelbrot(scheme)
    println("$(t), $(scheme), $(maximum(runtime_data.times))")
    update!(selector, maximum(runtime_data.times))
  end
end

main()
