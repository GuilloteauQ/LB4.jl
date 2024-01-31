include("../main.jl")

function diverges(c::Complex, nb_iterations::Int64, epsilon::Float64)::Bool
  z = 0 + im * 0
  i = 0
  while abs(z) < epsilon && i < nb_iterations
    z = z * z + c
    i = i + 1
  end
  i < nb_iterations
end

function mandelbrot()
  xmin = -4.0
  xmax =  4.0
  ymin = -1.0
  ymax =  1.0
  delta = 0.001


  nb_rows = floor(Int64, (ymax - ymin) / delta)
  nb_cols = floor(Int64, (xmax - xmin) / delta)

  data = Matrix{Bool}(undef, nb_rows, nb_cols)

  for x in 1:nb_cols
    for y in 1:nb_rows
      c = ((xmax - xmin) * x / nb_cols + xmin) + im * ((ymax - ymin) * y / nb_rows + ymin)
      data[y, x] = diverges(c, 10000, 100.0)
    end
  end
  data
end

function mandelbrot_par(scheme::Symbol)
  xmin = -4.0
  xmax =  4.0
  ymin = -1.0
  ymax =  1.0
  delta = 0.001


  nb_rows = floor(Int64, (ymax - ymin) / delta)
  nb_cols = floor(Int64, (xmax - xmin) / delta)

  data = Matrix{Bool}(undef, nb_rows, nb_cols)

  @lbthreads scheme for x in 1:nb_cols
    for y in 1:nb_rows
      c = ((xmax - xmin) * x / nb_cols + xmin) + im * ((ymax - ymin) * y / nb_rows + ymin)
      data[y, x] = diverges(c, 10000, 100.0)
    end
  end
  data
end

println("Mandelbrot --------------------")
schemes = [ :static, :gss, :tss, :dynamic, :fac2 ]
for scheme in schemes
  print(scheme)
  @time mandelbrot_par(scheme)
end
