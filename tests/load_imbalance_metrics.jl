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

function softmax(values::Array)
  out = exp.(values)
  s = sum(out)
  out / s
end

function normalize(values::Array)
  min_v = minimum(values)
  max_v = maximum(values)

  (values .- min_v) / (max_v - min_v)
end

function metrics(times::Array, nb_threads)
  digits = 3
  mean_time = sum(times) / nb_threads
  sd_time = sqrt(sum(map(x -> x^2, times .- mean_time)) / (nb_threads - 1))

  cov = round(sd_time / mean_time, digits=digits)

  max_time = maximum(times) 
  bound = round(1.96 * sd_time / sqrt(nb_threads), digits=digits)
  percentage_imbalance = round((nb_threads / (nb_threads-1)) * ((max_time - mean_time) / max_time) * 100, digits=digits)

  total_area = max_time * nb_threads 
  min_time = minimum(times)
  imbalance_area = sum(times .- min_time)
  area_metric = round(100 * imbalance_area / total_area, digits=digits)

  distance = round(100 * (max_time - sum(times) / nb_threads) / max_time, digits=digits)


  println("c.o.v.: $(cov), p.i.: $(percentage_imbalance)%, area: $(area_metric)%, distance: $(distance)%")
end

function mandelbrot_par(scheme::Symbol)
  xmin = -4.0
  xmax =  4.0
  ymin = -1.0
  ymax =  1.0
  delta = 0.002


  nb_rows = floor(Int64, (ymax - ymin) / delta)
  nb_cols = floor(Int64, (xmax - xmin) / delta)

  #data = Matrix{Bool}(undef, nb_rows, nb_cols)
  total_iters = nb_cols * nb_rows

  li = @LB4.lbthreads scheme for k in 0:(nb_cols * nb_rows - 2)
    y = div(k, nb_cols) 
    x = (k) % nb_cols
    c = ((xmax - xmin) * x / nb_cols + xmin) + im * ((ymax - ymin) * y / nb_rows + ymin)
    diverges(c, 10000, 2.0)
    #data[y + 1, x + 1] = diverges(c, 10000, 2.0)
  end
  println(li.times)
  metrics(li.times, li.nb_threads)
  soft_times = softmax(li.times)
  println(soft_times)
  metrics(soft_times, li.nb_threads)
  metrics(normalize(li.times), li.nb_threads)
end

mandelbrot_par(:fac2)

